import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import pandas as pd
import numpy as np

def threshold_zero(tensor, threshold=1e-3):
    """ Thresholds near-zero values to enforce sparsity. """
    return torch.where(tensor < threshold, torch.tensor(0.0, device=tensor.device), tensor)

# Modified VAE with presence decoder
class HybridVAE(nn.Module):
    def __init__(self, input_dim, latent_dim=32):
        super(HybridVAE, self).__init__()

        # Encoder
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.ReLU(True),
            nn.Linear(256, 128),
            nn.ReLU(True),
            nn.Linear(128, 64),
            nn.ReLU(True)
        )

        self.z_mean = nn.Linear(64, latent_dim)
        self.z_log_var = nn.Linear(64, latent_dim)

        # Decoder for abundance
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 64),
            nn.ReLU(True),
            nn.Linear(64, 128),
            nn.ReLU(True),
            nn.Linear(128, 256),
            nn.ReLU(True),
            nn.Linear(256, input_dim)
        )

        # Decoder for presence (Bernoulli probabilities)
        self.presence_decoder = nn.Sequential(
            nn.Linear(latent_dim, 64),
            nn.ReLU(True),
            nn.Linear(64, input_dim),
            nn.Sigmoid()
        )

    def encode(self, x):
        h = self.encoder(x)
        mu = self.z_mean(h)
        logvar = self.z_log_var(h)
        return mu, logvar

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z):
        return self.decoder(z)

    def decode_presence(self, z):
        return self.presence_decoder(z)

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon_x = self.decode(z)
        pres_prob = self.decode_presence(z)
        return recon_x, pres_prob, mu, logvar

def vae_loss(x, recon_x, pres_prob, mu, logvar,
                           lambda_pres, lambda_kl, lambda_abun,
                           gamma, lambda_moments, delta,
                           matching_param, SigmaHat, Mask, Sigma_norm,
                           target_corr):
    mask = (x > 0).float()
    batch_size = x.size(0)

    # VAE losses
    bernoulli = torch.distributions.Bernoulli(probs=pres_prob)
    log_prob_presence = bernoulli.log_prob(mask)
    loss_presence = -log_prob_presence.sum()

    loss_abundance = F.mse_loss(recon_x * mask, x, reduction='sum')
    kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())

    # Deep Knockoff: MMD
    n = x.shape[0] // 2
    X1, X2 = x[:n], x[n:2*n]
    Xk1, Xk2 = recon_x[:n], recon_x[n:2*n]

    Z1 = torch.cat((X1, Xk1), dim=1)
    Z2 = torch.cat((X2, Xk2), dim=1)
    Z3 = Z2.clone()

    p = x.shape[1]
    swap_inds = np.where(np.random.binomial(1, 0.5, size=p))[0]
    Z3[:, swap_inds] = Xk2[:, swap_inds]
    Z3[:, swap_inds + p] = X2[:, swap_inds]

    mmd_full = mix_rbf_mmd2_loss(Z1, Z2, matching_param)
    mmd_swap = mix_rbf_mmd2_loss(Z1, Z3, matching_param)

    # Deep Knockoff: second-order
    mX = x - x.mean(0, keepdim=True)
    mXk = recon_x - recon_x.mean(0, keepdim=True)
    SXkXk = torch.mm(mXk.T, mXk) / mXk.shape[0]
    SXXk = torch.mm(mX.T, mXk) / mXk.shape[0]
    loss_1m = (x.mean(0) - recon_x.mean(0)).pow(2).sum()
    loss_2m = ((SigmaHat - SXkXk).pow(2).sum() + (Mask * (SigmaHat - SXXk)).pow(2).sum()) / (SigmaHat.pow(2).sum() if Sigma_norm else 1.0)
    loss_moments = loss_1m + loss_2m

    # Deep Knockoff: decorrelation
    eps = 1e-3
    scaleX = mX.pow(2).mean(0, keepdim=True)
    scaleXk = mXk.pow(2).mean(0, keepdim=True)
    mXs = mX / (eps + torch.sqrt(scaleX))
    mXks = mXk / (eps + torch.sqrt(scaleXk))
    corr_XXk = (mXs * mXks).mean(0)
    loss_corr = (corr_XXk - target_corr).pow(2).mean()

    # Final loss
    total_loss = (
        lambda_pres * loss_presence +
        lambda_abun * loss_abundance +  # added lambda_abun
        lambda_kl * kl_loss +
        gamma * (mmd_full + mmd_swap) +
        lambda_moments * loss_moments +
        delta * loss_corr
    ) / batch_size

    return total_loss


def train_vae(
    model, otu_matrix,
    lambda_pres, lambda_kl, lambda_abun,
    gamma, lambda_moments, delta,
    matching_param, SigmaHat, Mask, target_corr,
    num_epochs=50000, batch_size=32,
    lr=0.01, weight_decay=1e-4,
    early_stopping_patience=200
):
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=10)

    dataset = TensorDataset(otu_matrix)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    best_loss = float('inf')
    epochs_no_improve = 0
    best_model_path = "best_vae_model.pth"

    for epoch in range(num_epochs):
        total_loss = 0
        for batch in dataloader:
            batch_data = batch[0]
            optimizer.zero_grad()
            recon_x, pres_prob, mu, logvar = model(batch_data)
            loss = vae_loss(batch_data, recon_x, pres_prob, mu, logvar,
                                   lambda_pres, lambda_kl, lambda_abun,
                                   gamma, lambda_moments, delta,
                                   matching_param, SigmaHat, Mask, True,
                                   target_corr)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

        avg_loss = total_loss / len(dataloader)
        scheduler.step(avg_loss)

        if avg_loss < best_loss:
            best_loss = avg_loss
            epochs_no_improve = 0
            torch.save(model.state_dict(), best_model_path)
        else:
            epochs_no_improve += 1

        if epoch % 200 == 0:
            print(f"Epoch {epoch}: Average Loss = {avg_loss:.4f}")

        if epochs_no_improve >= early_stopping_patience:
            print(f"Early stopping at epoch {epoch}, best loss = {best_loss:.4f}")
            break

    model.load_state_dict(torch.load(best_model_path))
    model.eval()
    print(f"Training complete! Best loss: {best_loss:.4f}")
    print(f"Best model saved at: {best_model_path}")
    return model, best_loss

def generate_samples(model, num_samples, taxa_names, latent_dim=32, device="cpu"):
    model.eval()
    model.to(device)
    z = torch.randn(num_samples, latent_dim).to(device)
    with torch.no_grad():
        abundance = model.decode(z)
        presence_prob = model.decode_presence(z)
        presence_sample = torch.bernoulli(presence_prob)
        generated_samples = abundance * presence_sample

    generated_samples = threshold_zero(generated_samples)
    generated_samples_np = generated_samples.cpu().numpy()
    df_generated = pd.DataFrame(generated_samples_np, columns=taxa_names)
    return df_generated

def covariance_diff_biased(X, Xk, SigmaHat, Mask, scale=1.0):
    """ Second-order loss function, as described in deep knockoffs manuscript
    :param X: input data
    :param Xk: generated knockoffs
    :param SigmaHat: target covariance matrix
    :param Mask: masking the diagonal of Cov(X,Xk)
    :param scale: scaling the loss function
    :return: second-order loss function
    """

    # Center X,Xk
    mX  = X  - torch.mean(X,0,keepdim=True)
    mXk = Xk - torch.mean(Xk,0,keepdim=True)
    # Compute covariance matrices
    SXkXk = torch.mm(torch.t(mXk),mXk)/mXk.shape[0]
    SXXk  = torch.mm(torch.t(mX),mXk)/mXk.shape[0]

    # Compute loss
    T  = (SigmaHat-SXkXk).pow(2).sum() / scale
    T += (Mask*(SigmaHat-SXXk)).pow(2).sum() / scale
    return T

# new testing function for diagnosis

def norm(X, p=2):
    if(p==np.inf):
        return(torch.max(torch.abs(X)))
    else:
        return(torch.norm(X,p))

def compute_diagnostics(x, recon_x):
    """
    Evaluate diagnostics between input x and reconstruction recon_x.
    Computes differences in moments and correlations.
    """
    diagnostics = {}

    # Mean difference
    D_mean = x.mean(0) - recon_x.mean(0)
    diagnostics["Mean"] = (D_mean * D_mean).mean().item()

    # Center X, Xk
    mX = x - x.mean(0, keepdim=True)
    mXk = recon_x - recon_x.mean(0, keepdim=True)

    # Scale
    scaleX = (mX * mX).mean(0, keepdim=True)
    scaleXk = (mXk * mXk).mean(0, keepdim=True)
    scaleX[scaleX <= 0] = 1.0
    scaleXk[scaleXk <= 0] = 1.0

    # Correlation between X and Xk
    mXs = mX / torch.sqrt(scaleX)
    mXks = mXk / torch.sqrt(scaleXk)
    diagnostics["Corr-Diag"] = (mXs * mXks).mean(0).mean().item()

    # Covariances
    # Sigma = torch.mm(mXs.T, mXs) / mXs.shape[0]
    # Sigma_ko = torch.mm(mXks.T, mXks) / mXks.shape[0]
    # DK_2 = norm(Sigma_ko - Sigma) / norm(Sigma)
    # diagnostics["Corr-Full"] = DK_2.item()

    # Corr-Swap: exclude diagonal
    # p = x.shape[1]
    # Mask = torch.ones(p, p) - torch.eye(p)
    # SigIntra_est = torch.mm(mX.T, mXk) / mXk.shape[0]
    # DS_2 = norm(Mask * (SigIntra_est - Sigma)) / norm(Sigma)
    # diagnostics["Corr-Swap"] = DS_2.item()

    # NEW: Add MMD-Full and MMD-Swap
    n = x.shape[0] // 2
    X1, X2 = x[:n], x[n:2*n]
    Xk1, Xk2 = recon_x[:n], recon_x[n:2*n]

    Z1 = torch.cat((X1, Xk1), dim=1)
    Z2 = torch.cat((X2, Xk2), dim=1)
    Z3 = Z2.clone()

    p = x.shape[1]
    swap_inds = np.where(np.random.binomial(1, 0.5, size=p))[0]
    Z3[:, swap_inds] = Xk2[:, swap_inds]
    Z3[:, swap_inds + p] = X2[:, swap_inds]

    matching_param = [1., 2., 4., 8., 16., 32., 64., 128.]
    diagnostics["MMD-Full"] = mix_rbf_mmd2_loss(Z1, Z2, matching_param).item()
    diagnostics["MMD-Swap"] = mix_rbf_mmd2_loss(Z1, Z3, matching_param).item()


    return diagnostics

min_var_est = 1e-8

# Consider linear time MMD with a linear kernel:
# K(f(x), f(y)) = f(x)^Tf(y)
# h(z_i, z_j) = k(x_i, x_j) + k(y_i, y_j) - k(x_i, y_j) - k(x_j, y_i)
#             = [f(x_i) - f(y_i)]^T[f(x_j) - f(y_j)]
#
# f_of_X: batch_size * k
# f_of_Y: batch_size * k
def linear_mmd2(f_of_X, f_of_Y):
    loss = 0.0
    delta = f_of_X - f_of_Y
    loss = torch.mean((delta[:-1] * delta[1:]).sum(1))
    return loss


# Consider linear time MMD with a polynomial kernel:
# K(f(x), f(y)) = (alpha*f(x)^Tf(y) + c)^d
# f_of_X: batch_size * k
# f_of_Y: batch_size * k
def poly_mmd2(f_of_X, f_of_Y, d=2, alpha=1.0, c=2.0):
    K_XX = (alpha * (f_of_X[:-1] * f_of_X[1:]).sum(1) + c)
    K_XX_mean = torch.mean(K_XX.pow(d))

    K_YY = (alpha * (f_of_Y[:-1] * f_of_Y[1:]).sum(1) + c)
    K_YY_mean = torch.mean(K_YY.pow(d))

    K_XY = (alpha * (f_of_X[:-1] * f_of_Y[1:]).sum(1) + c)
    K_XY_mean = torch.mean(K_XY.pow(d))

    K_YX = (alpha * (f_of_Y[:-1] * f_of_X[1:]).sum(1) + c)
    K_YX_mean = torch.mean(K_YX.pow(d))

    return K_XX_mean + K_YY_mean - K_XY_mean - K_YX_mean


def _mix_rbf_kernel(X, Y, sigma_list):
    assert(X.size(0) == Y.size(0))
    m = X.size(0)

    Z = torch.cat((X, Y), 0)
    ZZT = torch.mm(Z, Z.t())
    diag_ZZT = torch.diag(ZZT).unsqueeze(1)
    Z_norm_sqr = diag_ZZT.expand_as(ZZT)
    exponent = Z_norm_sqr - 2 * ZZT + Z_norm_sqr.t()

    K = 0.0
    for sigma in sigma_list:
        gamma = 1.0 / (2 * sigma**2)
        K += torch.exp(-gamma * exponent)

    return K[:m, :m], K[:m, m:], K[m:, m:], len(sigma_list)

def _mix_imq_kernel(X,
               Y,
               sigma_list):

    assert(X.size(0) == Y.size(0))
    m = X.size(0)
    h_dim = X.size(1)

    Z = torch.cat((X, Y), 0)
    ZZT = torch.mm(Z, Z.t())
    diag_ZZT = torch.diag(ZZT).unsqueeze(1)
    Z_norm_sqr = diag_ZZT.expand_as(ZZT)

    exponent = Z_norm_sqr - 2 * ZZT + Z_norm_sqr.t()

    K = 0.0
    for sigma in sigma_list:
        gamma = 2 * h_dim * 1.0 * sigma**2
        K += gamma / (gamma + exponent)

    return K[:m, :m], K[:m, m:], K[m:, m:], len(sigma_list)

def mix_imq_mmd2(X, Y, sigma_list, biased=True):
    K_XX, K_XY, K_YY, d = _mix_imq_kernel(X, Y, sigma_list)
    # return _mmd2(K_XX, K_XY, K_YY, const_diagonal=d, biased=biased)
    return _mmd2(K_XX, K_XY, K_YY, const_diagonal=False, biased=biased)

def mix_rbf_mmd2(X, Y, sigma_list, biased=True):
    K_XX, K_XY, K_YY, d = _mix_rbf_kernel(X, Y, sigma_list)
    # return _mmd2(K_XX, K_XY, K_YY, const_diagonal=d, biased=biased)
    return _mmd2(K_XX, K_XY, K_YY, const_diagonal=False, biased=biased)

def mix_rbf_mmd2_loss(X, Y, sigma_list, biased=True):

    mmd_dist_ref = mix_rbf_mmd2(X, Y, sigma_list, biased=True)
    return torch.sqrt(F.relu(mmd_dist_ref))

def mix_rbf_mmd2_unbiased_loss(X, Y, sigma_list):

    K_XX, K_XY, K_YY, d = _mix_rbf_kernel(X, Y, sigma_list)
    mmd_dist_ref = _mmd2_ignore_diagonals(K_XX, K_XY, K_YY, const_diagonal=False, biased=False)
    return torch.sqrt(F.relu(torch.abs(mmd_dist_ref)))

def mix_rbf_mmd2_and_ratio(X, Y, sigma_list, biased=True):
    K_XX, K_XY, K_YY, d = _mix_rbf_kernel(X, Y, sigma_list)
    # return _mmd2_and_ratio(K_XX, K_XY, K_YY, const_diagonal=d, biased=biased)
    return _mmd2_and_ratio(K_XX, K_XY, K_YY, const_diagonal=False, biased=biased)


################################################################################
# Helper functions to compute variances based on kernel matrices
################################################################################


def _mmd2_ignore_diagonals(K_XX, K_XY, K_YY, const_diagonal=False, biased=False):
    m = K_XX.size(0)    # assume X, Y are same shape

    # Get the various sums of kernels that we'll use
    # Kts drop the diagonal, but we don't need to compute them explicitly
    if const_diagonal is not False:
        diag_X = diag_Y = diag_XY = const_diagonal
        sum_diag_X = sum_diag_Y = sum_diag_XY = m * const_diagonal
    else:
        diag_X = torch.diag(K_XX)            # (m,)
        diag_Y = torch.diag(K_YY)            # (m,)
        diag_XY = torch.diag(K_XY)           # (m,)
        sum_diag_X = torch.sum(diag_X)
        sum_diag_Y = torch.sum(diag_Y)
        sum_diag_XY = torch.sum(diag_XY)


    Kt_XX_sums = K_XX.sum(dim=1) - diag_X             # \tilde{K}_XX * e = K_XX * e - diag_X
    Kt_YY_sums = K_YY.sum(dim=1) - diag_Y             # \tilde{K}_YY * e = K_YY * e - diag_Y
    K_XY_sums_0 = K_XY.sum(dim=0) - diag_XY                     # K_{XY}^T * e

    Kt_XX_sum = Kt_XX_sums.sum()                       # e^T * \tilde{K}_XX * e
    Kt_YY_sum = Kt_YY_sums.sum()                       # e^T * \tilde{K}_YY * e
    K_XY_sum = K_XY_sums_0.sum()                       # e^T * K_{XY} * e

    if biased:
        mmd2 = ((Kt_XX_sum + sum_diag_X) / (m * m)
            + (Kt_YY_sum + sum_diag_Y) / (m * m)
            - 2.0 * (K_XY_sum + sum_diag_XY) / (m * m))
    else:
        mmd2 = ((Kt_XX_sum ) / (m * (m-1))
            + (Kt_YY_sum) / (m * (m-1))
            - 2.0 * (K_XY_sum) / (m * (m-1)))

    return mmd2


def _mmd2(K_XX, K_XY, K_YY, const_diagonal=False, biased=False):
    m = K_XX.size(0)    # assume X, Y are same shape

    # Get the various sums of kernels that we'll use
    # Kts drop the diagonal, but we don't need to compute them explicitly
    if const_diagonal is not False:
        diag_X = diag_Y = const_diagonal
        sum_diag_X = sum_diag_Y = m * const_diagonal
    else:
        diag_X = torch.diag(K_XX)                       # (m,)
        diag_Y = torch.diag(K_YY)                       # (m,)
        sum_diag_X = torch.sum(diag_X)
        sum_diag_Y = torch.sum(diag_Y)

    Kt_XX_sums = K_XX.sum(dim=1) - diag_X             # \tilde{K}_XX * e = K_XX * e - diag_X
    Kt_YY_sums = K_YY.sum(dim=1) - diag_Y             # \tilde{K}_YY * e = K_YY * e - diag_Y
    K_XY_sums_0 = K_XY.sum(dim=0)                     # K_{XY}^T * e

    Kt_XX_sum = Kt_XX_sums.sum()                       # e^T * \tilde{K}_XX * e
    Kt_YY_sum = Kt_YY_sums.sum()                       # e^T * \tilde{K}_YY * e
    K_XY_sum = K_XY_sums_0.sum()                       # e^T * K_{XY} * e

    if biased:
        mmd2 = ((Kt_XX_sum + sum_diag_X) / (m * m)
            + (Kt_YY_sum + sum_diag_Y) / (m * m)
            - 2.0 * K_XY_sum / (m * m))
    else:
        mmd2 = (Kt_XX_sum / (m * (m - 1))
            + Kt_YY_sum / (m * (m - 1))
            - 2.0 * K_XY_sum / (m * m))

    return mmd2


def _mmd2_and_ratio(K_XX, K_XY, K_YY, const_diagonal=False, biased=False):
    mmd2, var_est = _mmd2_and_variance(K_XX, K_XY, K_YY, const_diagonal=const_diagonal, biased=biased)
    loss = mmd2 / torch.sqrt(torch.clamp(var_est, min=min_var_est))
    return loss, mmd2, var_est


def _mmd2_and_variance(K_XX, K_XY, K_YY, const_diagonal=False, biased=False):
    m = K_XX.size(0)    # assume X, Y are same shape

    # Get the various sums of kernels that we'll use
    # Kts drop the diagonal, but we don't need to compute them explicitly
    if const_diagonal is not False:
        diag_X = diag_Y = const_diagonal
        sum_diag_X = sum_diag_Y = m * const_diagonal
        sum_diag2_X = sum_diag2_Y = m * const_diagonal**2
    else:
        diag_X = torch.diag(K_XX)                       # (m,)
        diag_Y = torch.diag(K_YY)                       # (m,)
        sum_diag_X = torch.sum(diag_X)
        sum_diag_Y = torch.sum(diag_Y)
        sum_diag2_X = diag_X.dot(diag_X)
        sum_diag2_Y = diag_Y.dot(diag_Y)

    Kt_XX_sums = K_XX.sum(dim=1) - diag_X             # \tilde{K}_XX * e = K_XX * e - diag_X
    Kt_YY_sums = K_YY.sum(dim=1) - diag_Y             # \tilde{K}_YY * e = K_YY * e - diag_Y
    K_XY_sums_0 = K_XY.sum(dim=0)                     # K_{XY}^T * e
    K_XY_sums_1 = K_XY.sum(dim=1)                     # K_{XY} * e

    Kt_XX_sum = Kt_XX_sums.sum()                       # e^T * \tilde{K}_XX * e
    Kt_YY_sum = Kt_YY_sums.sum()                       # e^T * \tilde{K}_YY * e
    K_XY_sum = K_XY_sums_0.sum()                       # e^T * K_{XY} * e

    Kt_XX_2_sum = (K_XX ** 2).sum() - sum_diag2_X      # \| \tilde{K}_XX \|_F^2
    Kt_YY_2_sum = (K_YY ** 2).sum() - sum_diag2_Y      # \| \tilde{K}_YY \|_F^2
    K_XY_2_sum  = (K_XY ** 2).sum()                    # \| K_{XY} \|_F^2

    if biased:
        mmd2 = ((Kt_XX_sum + sum_diag_X) / (m * m)
            + (Kt_YY_sum + sum_diag_Y) / (m * m)
            - 2.0 * K_XY_sum / (m * m))
    else:
        mmd2 = (Kt_XX_sum / (m * (m - 1))
            + Kt_YY_sum / (m * (m - 1))
            - 2.0 * K_XY_sum / (m * m))

    var_est = (
        2.0 / (m**2 * (m - 1.0)**2) * (2 * Kt_XX_sums.dot(Kt_XX_sums) - Kt_XX_2_sum + 2 * Kt_YY_sums.dot(Kt_YY_sums) - Kt_YY_2_sum)
        - (4.0*m - 6.0) / (m**3 * (m - 1.0)**3) * (Kt_XX_sum**2 + Kt_YY_sum**2)
        + 4.0*(m - 2.0) / (m**3 * (m - 1.0)**2) * (K_XY_sums_1.dot(K_XY_sums_1) + K_XY_sums_0.dot(K_XY_sums_0))
        - 4.0*(m - 3.0) / (m**3 * (m - 1.0)**2) * (K_XY_2_sum) - (8 * m - 12) / (m**5 * (m - 1)) * K_XY_sum**2
        + 8.0 / (m**3 * (m - 1.0)) * (
            1.0 / m * (Kt_XX_sum + Kt_YY_sum) * K_XY_sum
            - Kt_XX_sums.dot(K_XY_sums_1)
            - Kt_YY_sums.dot(K_XY_sums_0))
        )
    return mmd2, var_est

