import numpy as np
import matplotlib.pyplot as plt

# Simulation parameters
n_dyads = 100
np.random.seed(42)

# Generate solo performances (S1, S2) from 0 to 1
S1 = np.random.beta(2, 2, n_dyads)  # Beta distribution for more variation in middle range
S2 = np.random.beta(2, 2, n_dyads)

# Simulate dyadic performances (D1, D2)
# Model: dyadic performance shows regression to mean + noise
# Better performers tend to drop slightly, worse performers tend to improve
# This creates convergence tendency with variability

def simulate_dyadic_performance(S1, S2):
    """
    Simulate dyadic performance with convergence/divergence dynamics
    """
    mean_performance = (S1 + S2) / 2
    
    # Social modulation: tendency toward partner's level + random variation
    convergence_strength = -0.3  # How much players converge (0=none, 1=complete)
    noise_level = 0.1
    
    D1 = S1 + convergence_strength * (mean_performance - S1) + np.random.normal(0, noise_level, len(S1))
    D2 = S2 + convergence_strength * (mean_performance - S2) + np.random.normal(0, noise_level, len(S2))
    
    # Clip to [0, 1] range
    D1 = np.clip(D1, 0, 1)
    D2 = np.clip(D2, 0, 1)
    
    return D1, D2

D1, D2 = simulate_dyadic_performance(S1, S2)

# Identify better and worse performers based on solo performance
Sbetter = np.maximum(S1, S2)
Sworse = np.minimum(S1, S2)

# Corresponding dyadic performances
Dbetter = np.where(S1 > S2, D1, D2)
Dworse = np.where(S1 > S2, D2, D1)

# Calculate metrics
solo_diff = Sbetter - Sworse  # Performance gap in solo condition
social_mod_diff = (Dbetter - Sbetter) - (Dworse - Sworse)  # Social modulation difference

# Plotting
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Main plot: Social modulation difference vs Solo difference
ax = axes[0, 0]
scatter = ax.scatter(solo_diff, social_mod_diff, alpha=0.3, s=20)
ax.axhline(y=0, color='k', linestyle='--', alpha=0.5, linewidth=1)
ax.set_xlabel('Solo Performance Gap (Sbetter - Sworse)', fontsize=11)
ax.set_ylabel('Social Modulation Difference\n(Dbetter - Sbetter) - (Dworse - Sworse)', fontsize=11)
ax.set_title('Social Modulation vs Solo Performance Gap', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)

# Add trend line
z = np.polyfit(solo_diff, social_mod_diff, 1)
p = np.poly1d(z)
x_trend = np.linspace(solo_diff.min(), solo_diff.max(), 100)
ax.plot(x_trend, p(x_trend), "r-", alpha=0.8, linewidth=2, label=f'Trend: y={z[0]:.3f}x+{z[1]:.3f}')

# Add theoretical perfect convergence line (D1 = D2)
# For perfect convergence: (Dbetter - Sbetter) - (Dworse - Sworse) = -(Sbetter - Sworse)
ax.plot(x_trend, -x_trend, "b-", alpha=0.8, linewidth=2, linestyle='--', label='Perfect convergence (D1=D2)')

# Add theoretical maximal divergence lines
# Upper bound: Dbetter=1, Dworse=0 → y = 1 - x
# Lower bound: Dbetter=0, Dworse=1 → y = -x - 1
ax.plot(x_trend, 1 - x_trend, "purple", alpha=0.6, linewidth=2, linestyle=':', label='Max divergence (Dbetter→1, Dworse→0)')
ax.plot(x_trend, -x_trend - 1, "orange", alpha=0.6, linewidth=2, linestyle=':', label='Max divergence (Dbetter→0, Dworse→1)')
ax.legend()

# Convergence vs Divergence histogram
ax = axes[0, 1]
convergence = (Dbetter - Dworse) < (Sbetter - Sworse)
conv_pct = np.mean(convergence) * 100
div_pct = 100 - conv_pct
ax.bar(['Convergence', 'Divergence'], [conv_pct, div_pct], color=['green', 'red'], alpha=0.7)
ax.set_ylabel('Percentage of Dyads (%)', fontsize=11)
ax.set_title('Convergence vs Divergence', fontsize=12, fontweight='bold')
ax.set_ylim([0, 100])
for i, (label, val) in enumerate(zip(['Convergence', 'Divergence'], [conv_pct, div_pct])):
    ax.text(i, val + 2, f'{val:.1f}%', ha='center', fontsize=10, fontweight='bold')

# Individual changes: better vs worse performers
ax = axes[1, 0]
better_change = Dbetter - Sbetter
worse_change = Dworse - Sworse
ax.scatter(better_change, worse_change, alpha=0.3, s=20)
ax.axhline(y=0, color='k', linestyle='--', alpha=0.5, linewidth=1)
ax.axvline(x=0, color='k', linestyle='--', alpha=0.5, linewidth=1)
ax.plot([-0.5, 0.5], [-0.5, 0.5], 'r--', alpha=0.5, linewidth=1, label='Equal change')
ax.set_xlabel('Better Performer Change (Dbetter - Sbetter)', fontsize=11)
ax.set_ylabel('Worse Performer Change (Dworse - Sworse)', fontsize=11)
ax.set_title('Individual Performance Changes', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()
ax.axis('equal')

# Distribution of social modulation difference
ax = axes[1, 1]
ax.hist(social_mod_diff, bins=50, alpha=0.7, edgecolor='black')
ax.axvline(x=0, color='r', linestyle='--', linewidth=2, label='No differential effect')
ax.axvline(x=np.mean(social_mod_diff), color='g', linestyle='-', linewidth=2, label=f'Mean = {np.mean(social_mod_diff):.3f}')
ax.set_xlabel('Social Modulation Difference', fontsize=11)
ax.set_ylabel('Frequency', fontsize=11)
ax.set_title('Distribution of Social Modulation Difference', fontsize=12, fontweight='bold')
ax.legend()

plt.tight_layout()
plt.savefig('cpr_convergence_analysis.png', dpi=300, bbox_inches='tight')
plt.show()

# Print summary statistics
print("="*60)
print("SIMULATION SUMMARY")
print("="*60)
print(f"Number of dyads: {n_dyads}")
print(f"\nConvergence: {conv_pct:.1f}%")
print(f"Divergence: {div_pct:.1f}%")
print(f"\nMean social modulation difference: {np.mean(social_mod_diff):.4f}")
print(f"Std social modulation difference: {np.std(social_mod_diff):.4f}")
print(f"\nCorrelation (solo gap vs social mod diff): {np.corrcoef(solo_diff, social_mod_diff)[0,1]:.4f}")
print(f"\nMean change for better performers: {np.mean(better_change):.4f}")
print(f"Mean change for worse performers: {np.mean(worse_change):.4f}")
print("="*60)

