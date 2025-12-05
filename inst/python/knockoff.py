import os
# Set threading environment variables BEFORE importing torch
# This prevents conflicts with R's threading model
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import pandas as pd
import numpy as np

# Configure PyTorch to use single thread to avoid conflicts with R
torch.set_num_threads(1)
torch.set_num_interop_threads(1)

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
        
        # Initialize weights to prevent extreme values that could cause NaN
        self._initialize_weights()
    
    def _initialize_weights(self):
        """Initialize weights to prevent numerical instability."""
        for m in self.modules():
            if isinstance(m, nn.Linear):
                # Use Xavier/Glorot initialization with smaller gain
                nn.init.xavier_uniform_(m.weight, gain=0.5)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0.0)

    def encode(self, x):
        # Check for NaN/Inf in input
        if torch.isnan(x).any() or torch.isinf(x).any():
            raise ValueError("Input x contains NaN or Inf values")
        
        h = self.encoder(x)
        
        # Check for NaN/Inf in encoder hidden state
        if torch.isnan(h).any() or torch.isinf(h).any():
            raise ValueError("Encoder output contains NaN or Inf values")
        
        mu = self.z_mean(h)
        logvar = self.z_log_var(h)
        
        # Clamp mu to prevent extreme values (though less critical than logvar)
        mu = torch.clamp(mu, min=-10.0, max=10.0)
        
        return mu, logvar

    def reparameterize(self, mu, logvar):
        # Clamp logvar to prevent extreme values that could cause NaN
        # logvar controls std = exp(0.5 * logvar), so clamping prevents extreme std
        logvar = torch.clamp(logvar, min=-10.0, max=10.0)
        std = torch.exp(0.5 * logvar)
        
        # Clamp std to prevent extreme values (safety check)
        std = torch.clamp(std, min=1e-6, max=1e3)
        
        eps = torch.randn_like(std)
        z = mu + eps * std
        
        # Clamp z to prevent extreme latent values
        z = torch.clamp(z, min=-50.0, max=50.0)
        
        return z

    def decode(self, z):
        return self.decoder(z)

    def decode_presence(self, z):
        # Check for NaN/Inf in latent z
        if torch.isnan(z).any() or torch.isinf(z).any():
            raise ValueError("Latent z contains NaN or Inf values")
        
        # Process through decoder layers manually to clamp logits before sigmoid
        # This prevents numerical overflow in sigmoid
        x = self.presence_decoder[0](z)  # Linear(latent_dim, 64)
        x = self.presence_decoder[1](x)  # ReLU
        logits = self.presence_decoder[2](x)  # Linear(64, input_dim)
        
        # Clamp logits before sigmoid to prevent numerical overflow
        # Sigmoid is stable for inputs in [-20, 20], but clamp more conservatively
        # This is the KEY fix: extreme logits cause sigmoid to produce NaN
        logits = torch.clamp(logits, min=-15.0, max=15.0)
        
        # Apply sigmoid (now safe from overflow)
        pres_prob = torch.sigmoid(logits)
        
        return pres_prob

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon_x = self.decode(z)
        pres_prob = self.decode_presence(z)
        
        # Clamp pres_prob to prevent NaN (sigmoid should already be in [0,1], but clamp for safety)
        eps = 1e-7
        pres_prob = torch.clamp(pres_prob, min=eps, max=1.0 - eps)
        
        # Check for NaN and replace if found
        if torch.isnan(pres_prob).any():
            pres_prob = torch.where(torch.isnan(pres_prob),
                                    torch.tensor(0.5, device=pres_prob.device, dtype=pres_prob.dtype),
                                    pres_prob)
            pres_prob = torch.clamp(pres_prob, min=eps, max=1.0 - eps)
        
        # Check for NaN in other outputs
        if torch.isnan(recon_x).any():
            recon_x = torch.where(torch.isnan(recon_x),
                                  torch.tensor(0.0, device=recon_x.device, dtype=recon_x.dtype),
                                  recon_x)
        
        # Return in format matching R: (recon_x, pres_prob, mu, logvar, input_x)
        return {
            'recon_x': recon_x,
            'pres_prob': pres_prob,
            'mu': mu,
            'logvar': logvar,
            'input_x': x
        }

def vae_loss(output, lambda_pres, lambda_kl, lambda_abun,
             gamma_full, gamma_swap, lambda_moments, delta_corr,
             sigma_list, SigmaHat, Mask, target_t):
    """
    Loss function matching R torch implementation exactly.
    """
    recon_x = output['recon_x']
    input_x = output['input_x']
    pres_prob = output['pres_prob']
    mu = output['mu']
    logvar = output['logvar']
    
    # Clamp pres_prob to valid range [eps, 1-eps] to avoid NaN and ensure valid probabilities
    # This prevents numerical instability in sigmoid and Bernoulli distribution
    eps = 1e-7
    pres_prob = torch.clamp(pres_prob, min=eps, max=1.0 - eps)
    
    # Check for NaN values and replace with default value if found
    if torch.isnan(pres_prob).any():
        # Replace NaN with 0.5 (neutral probability)
        pres_prob = torch.where(torch.isnan(pres_prob), 
                                torch.tensor(0.5, device=pres_prob.device, dtype=pres_prob.dtype),
                                pres_prob)
        # Re-clamp after NaN replacement
        pres_prob = torch.clamp(pres_prob, min=eps, max=1.0 - eps)
    
    mask = (input_x > 0).float()
    n_batch = input_x.size(0)
    
    # Presence loss (Bernoulli)
    bernoulli = torch.distributions.Bernoulli(probs=pres_prob)
    loss_pres = -bernoulli.log_prob(mask).sum()
    
    # Abundance loss: use xk = recon_x * pres_prob (not recon_x * mask)
    xk = recon_x * pres_prob
    loss_abun = F.mse_loss(xk, input_x, reduction='sum')
    
    # KL loss
    kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - torch.exp(logvar))
    
    # Deep Knockoff: MMD
    # Need at least 2 samples for MMD (1 for each group)
    if n_batch < 2:
        # If batch is too small, set MMD losses to zero
        mmd_full = torch.tensor(0.0, device=input_x.device, dtype=torch.float32)
        mmd_swap = torch.tensor(0.0, device=input_x.device, dtype=torch.float32)
    else:
        n_half = n_batch // 2
        # Ensure we have at least 1 sample in each half
        if n_half < 1:
            n_half = 1
        
        # Make sure we don't exceed bounds
        end_idx = min(2 * n_half, n_batch)
        X1 = input_x[:n_half]
        X2 = input_x[n_half:end_idx]
        Xk1 = xk[:n_half]
        Xk2 = xk[n_half:end_idx]
        
        # Ensure X1 and X2 have same number of samples for MMD
        min_samples = min(X1.size(0), X2.size(0))
        if min_samples > 0:
            X1 = X1[:min_samples]
            X2 = X2[:min_samples]
            Xk1 = Xk1[:min_samples]
            Xk2 = Xk2[:min_samples]
            
            # Z1 = [X1, Xk1], Z2 = [Xk2, X2], Z3 = [X2, Xk2] with swaps
            Z1 = torch.cat([X1, Xk1], dim=1)
            Z2 = torch.cat([Xk2, X2], dim=1)
            Z3 = torch.cat([X2, Xk2], dim=1)
            
            p = input_x.size(1)
            swap = np.where(np.random.binomial(1, 0.5, size=p) == 1)[0]
            if len(swap) > 0:
                Z3[:, swap] = Xk2[:, swap]
                Z3[:, swap + p] = X2[:, swap]
            
            mmd_full = mix_rbf_mmd2_loss(Z1, Z2, sigma_list)
            mmd_swap = mix_rbf_mmd2_loss(Z1, Z3, sigma_list)
        else:
            mmd_full = torch.tensor(0.0, device=input_x.device, dtype=torch.float32)
            mmd_swap = torch.tensor(0.0, device=input_x.device, dtype=torch.float32)
    
    # Deep Knockoff: second-order moments
    # R uses mean(1, keepdim=True) which is mean along rows (dim=0 in PyTorch)
    mX = input_x - input_x.mean(0, keepdim=True)
    mXk = xk - xk.mean(0, keepdim=True)
    SXk = torch.mm(mXk.t(), mXk) / n_batch
    SXXk = torch.mm(mX.t(), mXk) / n_batch
    loss_1m = (input_x.mean(0) - xk.mean(0)).pow(2).sum()
    scale_f = torch.sum(SigmaHat.pow(2))
    loss_2m = ((SigmaHat - SXk).pow(2).sum() +
                (Mask * (SigmaHat - SXXk)).pow(2).sum()) / scale_f
    loss_mom = loss_1m + loss_2m
    
    # Deep Knockoff: decorrelation
    scaleX = mX.pow(2).mean(0, keepdim=True)
    scaleXk = mXk.pow(2).mean(0, keepdim=True)
    corr = ((mX / (torch.sqrt(scaleX) + 1e-3)) *
            (mXk / (torch.sqrt(scaleXk) + 1e-3))).mean(0)
    loss_corr = (corr - target_t).pow(2).mean()
    
    # Final loss
    loss = (lambda_pres * loss_pres +
            lambda_abun * loss_abun +
            lambda_kl * kl_loss +
            gamma_full * mmd_full +
            gamma_swap * mmd_swap +
            lambda_moments * loss_mom +
            delta_corr * loss_corr) / n_batch
    
    return loss


def train_vae(
    model, x,
    lambda_pres, lambda_kl, lambda_abun,
    gamma_full, gamma_swap, lambda_moments, delta_corr,
    sigma_list, SigmaHat, Mask, target_t,
    epochs=100, batch_size=50,
    lr=1e-3, weight_decay=1e-2,
    progress=True
):
    """
    Training function matching R torch implementation.
    Uses AdamW optimizer and trains for specified epochs.
    """
    optimizer = optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    
    dataset = TensorDataset(x)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    model.train()
    for epoch in range(epochs):
        total_loss = 0
        for batch in dataloader:
            batch_data = batch[0]
            optimizer.zero_grad()
            output = model(batch_data)
            
            # Check for NaN in output before computing loss
            if torch.isnan(output['pres_prob']).any() or torch.isnan(output['recon_x']).any():
                # Skip this batch if NaN detected (likely numerical instability)
                continue
            
            loss = vae_loss(output, lambda_pres, lambda_kl, lambda_abun,
                           gamma_full, gamma_swap, lambda_moments, delta_corr,
                           sigma_list, SigmaHat, Mask, target_t)
            
            # Check for NaN or Inf in loss
            if torch.isnan(loss) or torch.isinf(loss):
                # Skip this batch if loss is invalid
                continue
            
            loss.backward()
            
            # Gradient clipping to prevent exploding gradients
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            
            optimizer.step()
            total_loss += loss.item()

        if progress and (epoch + 1) % 10 == 0:
            avg_loss = total_loss / len(dataloader)
            print(f"Epoch {epoch + 1}/{epochs}: Average Loss = {avg_loss:.6f}")

    return model

def VAE_func_DK(x,
                 latent_dim=32,
                 lambda_kl=0,
                 lambda_abun=0,
                 lambda_pres=0,
                 gamma_full=1,
                 gamma_swap=1,
                 lambda_moments=0,
                 delta_corr=0,
                 sigma_list=[1, 2, 4, 8, 16, 32, 64, 128],
                 epochs=100,
                 batch_size=50,
                 lr=1e-3,
                 weight_decay=1e-2,
                 seed=123,
                 progress=True):
    """
    Main VAE function matching R torch VAE_func_DK implementation.
    
    Parameters:
    -----------
    x : numpy array or torch tensor
        Input data matrix (n_samples, n_features)
    latent_dim : int
        Dimension of latent space
    lambda_kl : float
        Weight for KL divergence loss
    lambda_abun : float
        Weight for abundance loss
    lambda_pres : float
        Weight for presence loss
    gamma_full : float
        Weight for MMD full loss
    gamma_swap : float
        Weight for MMD swap loss
    lambda_moments : float
        Weight for moment matching loss
    delta_corr : float
        Weight for correlation loss
    sigma_list : list
        List of sigma values for MMD kernel
    epochs : int
        Number of training epochs
    batch_size : int
        Batch size for training
    lr : float
        Learning rate
    weight_decay : float
        Weight decay for optimizer
    seed : int
        Random seed
    
    Returns:
    --------
    dict with keys:
        - knockoff_x : numpy array of generated knockoffs
        - recon_x : numpy array of reconstructions
    """
    # Set random seeds
    torch.manual_seed(seed)
    np.random.seed(seed)
    
    # Convert input to torch tensor if needed
    if isinstance(x, np.ndarray):
        x_tensor = torch.tensor(x, dtype=torch.float32)
    else:
        x_tensor = x.float()
    
    p = x_tensor.shape[1]
    
    # Compute covariance matrix and mask (matching R code)
    # R: cov(x) computes covariance of columns (features), giving p×p matrix
    # NumPy: np.cov(x.T) does the same for (n_samples, n_features) input
    x_np = x_tensor.numpy() if isinstance(x_tensor, torch.Tensor) else x_tensor
    SigmaHat_t = torch.tensor(np.cov(x_np.T), dtype=torch.float32)
    Mask_t = torch.tensor(np.ones((p, p)) - np.eye(p), dtype=torch.float32)
    target_t = torch.zeros(p, dtype=torch.float32)
    
    # Ensure batch_size doesn't exceed dataset size
    n_samples = x_tensor.shape[0]
    effective_batch_size = min(batch_size, n_samples)
    if effective_batch_size < batch_size:
        print(f"Warning: batch_size ({batch_size}) exceeds dataset size ({n_samples}). Using batch_size={effective_batch_size}")
    
    # Initialize model
    model = HybridVAE(input_dim=p, latent_dim=latent_dim)
    
    # Train model
    model = train_vae(
        model, x_tensor,
        lambda_pres, lambda_kl, lambda_abun,
        gamma_full, gamma_swap, lambda_moments, delta_corr,
        sigma_list, SigmaHat_t, Mask_t, target_t,
        epochs=epochs, batch_size=effective_batch_size,
        lr=lr, weight_decay=weight_decay,
        progress=progress
    )
    
    # Generate knockoffs (matching R code)
    model.eval()
    with torch.no_grad():
        output = model(x_tensor)
    
    # Convert to numpy arrays (matching R as_array)
    input_x = output['input_x'].numpy()
    recon_x = output['recon_x'].numpy()
    pres_prob = output['pres_prob'].numpy()
    
    # Sample from Bernoulli distribution (matching R torch_bernoulli)
    mask = np.random.binomial(1, pres_prob).astype(float)
    
    # Generate knockoffs (matching R: knockoff_x <- recon_x * mask)
    knockoff_x = recon_x * mask
    
    return {
        'knockoff_x': knockoff_x,
        'recon_x': recon_x
    }

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

def mix_rbf_mmd2_loss(X, Y, sigma_list):
    """
    MMD loss using pairwise distances (matching R torch implementation).
    Uses torch.cdist to compute pairwise distances, then applies Gaussian kernels.
    """
    # Compute pairwise squared distances (p=2 norm, squared)
    XX = torch.cdist(X, X, p=2).pow(2)  # 2-norm distances ^ 2
    YY = torch.cdist(Y, Y, p=2).pow(2)
    XY = torch.cdist(X, Y, p=2).pow(2)
    
    mmd2 = torch.tensor(0.0, dtype=torch.float32, device=X.device)
    
    for s in sigma_list:
        g = 1.0 / (2 * s**2)  # bandwidth parameter
        mmd2 = mmd2 + (
            torch.mean(torch.exp(-g * XX)) +  # Gaussian kernel
            torch.mean(torch.exp(-g * YY)) -
            2 * torch.mean(torch.exp(-g * XY))
        )
    
    return torch.sqrt(F.relu(mmd2))

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

