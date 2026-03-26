"""
Numerical Quantum Mechanics: Basis Expansion Methods
=====================================================
Lecture Example Code

Part 1: Radial Harmonic Oscillator in Coordinate Space (Matrix Eigenvalue Method)
Part 2: 1-Body Bound State Problem using Harmonic Oscillator Basis Expansion

Physics conventions (natural units): hbar = m = omega_HO = 1

Author: Graduate Lecture in Numerical Quantum Mechanics
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.special import factorial, genlaguerre, gamma
import matplotlib.gridspec as gridspec

# ============================================================
# PART 1: RADIAL HARMONIC OSCILLATOR IN COORDINATE SPACE
# ============================================================
#
# Radial Schrodinger equation (reduced wavefunction u = r*R):
#
#   [-1/2  d^2/dr^2  +  l(l+1)/(2r^2)  +  V(r)] u(r)  =  E u(r)
#
# For harmonic oscillator:  V(r) = 1/2 * r^2
#
# Boundary conditions:  u(0) = 0,  u(r_max) = 0  (Dirichlet)
#
# Finite difference for d^2/dr^2 (3-point central difference):
#   d^2u/dr^2 at point i  ~  [u_{i-1} - 2u_i + u_{i+1}] / dr^2
#
# Exact eigenvalues:  E_{nl} = 2n + l + 3/2,  n = 0, 1, 2, ...
# ============================================================


def build_hamiltonian_fd(r, V_func, l=0):
    """
    Build the 1D radial Hamiltonian matrix using 3-point finite differences.

    Grid: uniform, r[i] = (i+1)*dr, with Dirichlet BCs u(0)=u(r_max)=0.
    The boundary points are NOT included in r.

    Parameters
    ----------
    r      : 1D array of interior radial grid points
    V_func : callable, V_func(r) returns potential values
    l      : orbital angular momentum quantum number

    Returns
    -------
    H  : (N x N) Hamiltonian matrix  (real symmetric)
    """
    N  = len(r)
    dr = r[1] - r[0]

    # Kinetic energy matrix: T = -1/2 * d^2/dr^2
    #   T_{ii}     = +1/dr^2          (diagonal)
    #   T_{i,i+1}  = T_{i,i-1} = -1/(2 dr^2)  (off-diagonal)
    T_diag    =  np.ones(N)   / dr**2
    T_offdiag = -np.ones(N-1) / (2.0 * dr**2)

    T = (np.diag(T_diag)
         + np.diag(T_offdiag, k=+1)
         + np.diag(T_offdiag, k=-1))

    # Potential: centrifugal barrier + external potential
    V_centrifugal = l * (l + 1) / (2.0 * r**2)
    V_ext         = V_func(r)
    V = np.diag(V_centrifugal + V_ext)

    return T + V


def exact_ho_energy(n, l):
    """
    Exact radial HO eigenvalue (hbar=m=omega=1).
    E_{nl} = 2n + l + 3/2,   n = 0, 1, 2, ...
    """
    return 2.0*n + l + 1.5


def exact_ho_wavefunction(r, n, l):
    """
    Exact normalized radial HO wavefunction u_{nl}(r) = r * R_{nl}(r).

    u_{nl}(r) = N_{nl} * r^{l+1} * exp(-r^2/2) * L_n^{l+1/2}(r^2)

    where L_n^alpha is the associated Laguerre polynomial and
    N_{nl} is the normalization constant satisfying int |u|^2 dr = 1.
    """
    alpha = l + 0.5
    norm  = np.sqrt(2.0 * factorial(n) / gamma(n + alpha + 1))
    L     = genlaguerre(n, alpha)
    u     = norm * r**(l+1) * np.exp(-r**2 / 2.0) * L(r**2)
    return u


def V_harmonic_oscillator(r, omega=1.0):
    """3D isotropic HO potential: V(r) = 1/2 * omega^2 * r^2"""
    return 0.5 * omega**2 * r**2


def solve_ho_coordinate_space(l=0, r_max=8.0, N_grid=400, n_states=5):
    """
    Solve the radial HO eigenvalue problem in coordinate space
    by direct matrix diagonalization.

    Grid: r_i = i * dr  for i = 1, 2, ..., N_grid
    Dirichlet BCs:  u(r=0) = 0  (implicit),  u(r_max) = 0  (implicit).

    Parameters
    ----------
    l       : angular momentum quantum number
    r_max   : box radius (should be >> classical turning point)
    N_grid  : number of interior grid points
    n_states: number of states to return

    Returns
    -------
    r     : grid array
    E_num : numerical eigenvalues
    u_num : normalized eigenvectors (columns)
    """
    dr = r_max / (N_grid + 1)
    r  = np.arange(1, N_grid + 1) * dr    # interior: dr, 2dr, ..., N*dr

    H = build_hamiltonian_fd(r, V_harmonic_oscillator, l=l)
    E_num, u_num = eigh(H, subset_by_index=[0, n_states - 1])

    # Normalize: int_{0}^{r_max} |u(r)|^2 dr = 1
    for i in range(n_states):
        norm = np.sqrt(np.sum(u_num[:, i]**2) * dr)
        u_num[:, i] /= norm*np.sign(u_num[0,i])
        #if u_num[np.argmax(np.abs(u_num[:, i])), i] < 0:
        #    u_num[:, i] *= -1

    return r, E_num, u_num


# ============================================================
# PART 2: HARMONIC OSCILLATOR BASIS EXPANSION
# ============================================================
#
# Express the wavefunction as:
#       u(r) = sum_{n=0}^{N_basis-1}  c_n  * phi_n(r)
#
# where phi_n(r) = u_{nl}(r) are exact HO eigenfunctions (orthonormal).
#
# Hamiltonian in HO basis:
#       H_{ij} = <phi_i | T + V_centrifugal + V_ext | phi_j>
#
# Since phi_n are eigenstates of H_HO = T + V_centrifugal + V_HO:
#       <phi_i | T + V_centrifugal | phi_j> = E_i^HO * delta_{ij}
#                                             - <phi_i | V_HO | phi_j>
#
# Therefore:
#       H_{ij} = E_i^HO * delta_{ij}  +  <phi_i | V_ext - V_HO | phi_j>
#
# The secular equation  H * c = E * c  gives bound state energies.
# ============================================================


def compute_potential_matrix(r, basis_funcs, V_func):
    """
    Compute potential matrix elements  V_{ij} = <phi_i | V | phi_j>
    by numerical integration (trapezoidal rule).
    """
    N_basis = len(basis_funcs)
    V_vals  = V_func(r)
    V_mat   = np.zeros((N_basis, N_basis))

    for i in range(N_basis):
        for j in range(i, N_basis):
            Vij = np.trapezoid(basis_funcs[i] * V_vals * basis_funcs[j], r)
            V_mat[i, j] = Vij
            V_mat[j, i] = Vij

    return V_mat


def solve_with_ho_basis(V_ext_func, l=0, N_basis=20,
                        r_max=15.0, N_grid=2000, n_states=5):
    """
    Solve a 1-body radial Schrodinger equation using HO basis expansion.

    Steps:
      1. Construct HO basis functions on a fine grid.
      2. Build H_{ij} = E_i^HO*delta_{ij} + <phi_i|V_ext - V_HO|phi_j>
      3. Diagonalize H to get eigenvalues and expansion coefficients.

    Parameters
    ----------
    V_ext_func : callable, V_ext(r) — external potential
    l          : angular momentum
    N_basis    : number of HO basis functions
    r_max      : integration cutoff
    N_grid     : number of integration grid points
    n_states   : number of states to return

    Returns
    -------
    E_vals : lowest n_states eigenvalues
    E_vecs : corresponding eigenvectors (columns = expansion coefficients c_n)
    r      : integration grid
    """
    r = np.linspace(1e-6, r_max, N_grid)

    # Build orthonormal HO basis
    basis = np.array([exact_ho_wavefunction(r, n, l) for n in range(N_basis)])

    # Diagonal: HO kinetic + centrifugal eigenvalues
    HO_energies = np.array([exact_ho_energy(n, l) for n in range(N_basis)])

    # Off-diagonal: residual potential V_ext - V_HO
    V_ext_mat = compute_potential_matrix(r, basis, V_ext_func)
    V_ho_mat  = compute_potential_matrix(r, basis, V_harmonic_oscillator)

    H = np.diag(HO_energies) + (V_ext_mat - V_ho_mat)

    E_vals, E_vecs = eigh(H)

    return E_vals[:n_states], E_vecs[:, :n_states], r


def reconstruct_wavefunction(r, coeffs, l, N_basis):
    """
    Reconstruct u(r) = sum_{n} c_n * phi_n(r) from expansion coefficients.
    """
    u = np.zeros(len(r))
    for n in range(N_basis):
        u += coeffs[n] * exact_ho_wavefunction(r, n, l)
    if u[np.argmax(np.abs(u))] < 0:
        u *= -1
    return u


# ============================================================
# EXAMPLE POTENTIALS
# ============================================================

def V_woods_saxon(r, V0=50.0, R=5.0, a=0.5):
    """
    Woods-Saxon potential (nuclear physics):
        V(r) = -V0 / (1 + exp((r - R)/a))
    """
    return -V0 / (1.0 + np.exp((r - R) / a))


def V_morse(r, De=10.0, alpha=1.0, re=2.0):
    """
    Morse potential (molecular physics):
        V(r) = De * (1 - exp(-alpha*(r - re)))^2 - De
    """
    return De * (1.0 - np.exp(-alpha * (r - re)))**2 - De


# ============================================================
# PRINT SUMMARY TABLE
# ============================================================

def print_summary():
    print("=" * 70)
    print(" PART 1: HO Coordinate Space — Finite Difference (l=0, N=400, r_max=8)")
    print("=" * 70)
    print(f"{'n':>4}  {'E_numerical':>14}  {'E_exact':>10}  {'|Error|':>12}  {'Rel.Err.':>10}")
    print("-" * 58)
    r, E_num, _ = solve_ho_coordinate_space(l=0, r_max=8, N_grid=400, n_states=6)
    for n in range(6):
        E_ex = exact_ho_energy(n, l=0)
        err  = abs(E_num[n] - E_ex)
        print(f"{n:>4}  {E_num[n]:>14.8f}  {E_ex:>10.4f}  {err:>12.2e}  {err/E_ex:>10.2e}")

    print()
    print("=" * 70)
    print(" PART 2: HO Basis Expansion — Woods-Saxon Potential (l=0)")
    print("=" * 70)
    print(f"{'N_basis':>8}  {'E0':>12}  {'E1':>12}  {'E2':>12}")
    print("-" * 50)
    for N_b in [5, 8, 10, 12, 15, 20, 25, 30]:
        ev, _, _ = solve_with_ho_basis(V_woods_saxon, l=0, N_basis=N_b, n_states=3)
        print(f"{N_b:>8}  {ev[0]:>12.6f}  {ev[1]:>12.6f}  {ev[2]:>12.6f}")

    print()
    print("=" * 70)
    print(" PART 2: HO Basis Expansion — Morse Potential (bound states, l=0)")
    print("=" * 70)
    for N_b in [5, 10, 15, 20]:
        ev, _, _ = solve_with_ho_basis(V_morse, l=0, N_basis=N_b, n_states=6)
        bound = [f"{e:.5f}" for e in ev if e < 0]
        print(f"  N_basis={N_b:2d}: {bound}")


# ============================================================
# PLOT ALL RESULTS
# ============================================================

def plot_all():
    fig = plt.figure(figsize=(16, 20))
    fig.patch.set_facecolor('white')
    gs  = gridspec.GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.35)
    COLORS = ['#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#ff7f0e']

    # ── Panel 1A: HO wavefunctions in coordinate space ───────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    r, E_num, u_num = solve_ho_coordinate_space(l=0, r_max=8, N_grid=500, n_states=5)
    r_exact = np.linspace(0.1, 7.5, 80)
    for i in range(4):
        E_ex = exact_ho_energy(i, l=0)
        u_ex = exact_ho_wavefunction(r_exact, i, l=0)
        ax1.plot(r_exact, u_ex, '*', color=COLORS[i], lw=1.2, alpha=0.7)
        ax1.plot(r, u_num[:, i], '-',  color=COLORS[i], lw=2.2,
                 label=f'$n={i}$: $E_{{\\rm num}}={E_num[i]:.5f}$, $E_{{\\rm ex}}={E_ex:.4f}$')
    ax1.set_xlabel('$r$ (nat. units)', fontsize=12)
    ax1.set_ylabel('$u_{n,l=0}(r)$', fontsize=12)
    ax1.set_title('Part 1: HO Radial Wavefunctions ($l=0$)\n'
                  'Solid = FD matrix, Dashed = analytic', fontsize=11, fontweight='bold')
    ax1.legend(fontsize=8.5, loc='upper right'); ax1.set_xlim(0, 7.5)
    ax1.axhline(0, color='k', lw=0.5); ax1.grid(True, alpha=0.3)

    # ── Panel 1B: Eigenvalue error vs grid size (log-log) ────────────────
    ax2 = fig.add_subplot(gs[0, 1])
    grid_sizes = [50, 100, 200, 400, 800, 1600]
    for n_state in range(3):
        errs = []
        for N in grid_sizes:
            _, ev, _ = solve_ho_coordinate_space(l=0, r_max=8, N_grid=N, n_states=n_state+1)
            errs.append(abs(ev[n_state] - exact_ho_energy(n_state, 0)))
        ax2.loglog(grid_sizes, errs, 'o-', color=COLORS[n_state], lw=1.8, label=f'$n={n_state}$')
    h = 8.0 / np.array(grid_sizes)
    ax2.loglog(grid_sizes, 0.03 * h**2, 'k:', lw=1.5, label='$O(h^2)$')
    ax2.set_xlabel('Number of grid points $N$', fontsize=12)
    ax2.set_ylabel('$|E_{\\rm num} - E_{\\rm exact}|$', fontsize=12)
    ax2.set_title('Part 1: FD Eigenvalue Error vs Grid Size\n($\\propto h^2$ convergence expected)',
                  fontsize=11, fontweight='bold')
    ax2.legend(fontsize=10); ax2.grid(True, alpha=0.3, which='both')

    # ── Panel 2A: Basis convergence for Woods-Saxon ───────────────────────
    ax3 = fig.add_subplot(gs[1, 0])
    basis_sizes = [4, 6, 8, 10, 12, 15, 20, 25, 30]
    E_ref, _, _ = solve_with_ho_basis(V_woods_saxon, l=0, N_basis=40, n_states=3)
    for i_state, ls in enumerate(['-o', '-s', '-^']):
        errs = []
        for nb in basis_sizes:
            ev, _, _ = solve_with_ho_basis(V_woods_saxon, l=0, N_basis=nb, n_states=i_state+1)
            errs.append(abs(ev[i_state] - E_ref[i_state]))
        ax3.semilogy(basis_sizes, errs, ls, color=COLORS[i_state], lw=1.8, ms=7,
                     label=f'State $n={i_state}$')
    ax3.set_xlabel('Basis size $N_{\\rm basis}$', fontsize=12)
    ax3.set_ylabel('$|E - E_{\\rm ref}|$ (nat. units)', fontsize=12)
    ax3.set_title('Part 2: Basis Convergence\nWoods-Saxon ($l=0$)', fontsize=11, fontweight='bold')
    ax3.legend(fontsize=10); ax3.grid(True, alpha=0.3, which='both')
    ax3.text(0.52, 0.88, f'$E_0^{{\\rm ref}} = {E_ref[0]:.4f}$',
             transform=ax3.transAxes, fontsize=10,
             bbox=dict(boxstyle='round', fc='lightyellow'))

    # ── Panel 2B: Woods-Saxon wavefunctions ──────────────────────────────
    ax4 = fig.add_subplot(gs[1, 1])
    N_basis = 20
    E_ws, c_ws, r_ws = solve_with_ho_basis(V_woods_saxon, l=0, N_basis=N_basis,
                                            r_max=15, n_states=4)
    for i in range(3):
        u = reconstruct_wavefunction(r_ws, c_ws[:, i], l=0, N_basis=N_basis)
        ax4.plot(r_ws, u, color=COLORS[i], lw=1.8, label=f'$n={i}$,  $E={E_ws[i]:.3f}$')
    V_ws = V_woods_saxon(r_ws)
    ax4.plot(r_ws, V_ws / abs(V_ws.min()) * 0.7, 'k--', lw=1.4, alpha=0.6,
             label='$V_{\\rm WS}$ (rescaled)')
    ax4.set_xlabel('$r$', fontsize=12); ax4.set_ylabel('$u_{n,l=0}(r)$', fontsize=12)
    ax4.set_title(f'Part 2: Woods-Saxon Bound States\n(HO basis, $N_{{\\rm basis}}={N_basis}$, $l=0$)',
                  fontsize=11, fontweight='bold')
    ax4.legend(fontsize=9); ax4.set_xlim(0, 12)
    ax4.axhline(0, color='k', lw=0.5); ax4.grid(True, alpha=0.3)

    # ── Panel 3A: Morse potential wavefunctions ───────────────────────────
    ax5 = fig.add_subplot(gs[2, 0])
    E_morse, c_morse, r_morse = solve_with_ho_basis(V_morse, l=0, N_basis=20,
                                                     r_max=10, n_states=6)
    n_plot = 0
    for i in range(len(E_morse)):
        if E_morse[i] < 0 and n_plot < 4:
            u = reconstruct_wavefunction(r_morse, c_morse[:, i], l=0, N_basis=20)
            ax5.plot(r_morse, u, color=COLORS[n_plot], lw=1.8,
                     label=f'$n={n_plot}$, $E={E_morse[i]:.3f}$')
            n_plot += 1
    V_m = V_morse(r_morse)
    ax5.plot(r_morse, V_m / abs(V_m.min()) * 0.7, 'k--', lw=1.4, alpha=0.6,
             label='$V_{\\rm Morse}$ (rescaled)')
    ax5.set_xlabel('$r$ (nat. units)', fontsize=12); ax5.set_ylabel('$u_{n,l=0}(r)$', fontsize=12)
    ax5.set_title('Part 2: Morse Potential Bound States\n(HO basis, $l=0$)',
                  fontsize=11, fontweight='bold')
    ax5.legend(fontsize=9); ax5.set_xlim(0, 9)
    ax5.axhline(0, color='k', lw=0.5); ax5.grid(True, alpha=0.3)

    # ── Panel 3B: Expansion coefficients |c_n|^2 ─────────────────────────
    ax6 = fig.add_subplot(gs[2, 1])
    n_arr = np.arange(N_basis)
    for i in range(3):
        ax6.bar(n_arr + i*0.25, c_ws[:, i]**2, width=0.25,
                color=COLORS[i], alpha=0.80, label=f'State $n={i}$')
    ax6.set_xlabel('Basis index $n$', fontsize=12)
    ax6.set_ylabel('$|c_n|^2$', fontsize=12)
    ax6.set_title('Part 2: HO Expansion Coefficients $|c_n|^2$\nWoods-Saxon Bound States',
                  fontsize=11, fontweight='bold')
    ax6.legend(fontsize=10); ax6.grid(True, alpha=0.3, axis='y')

    fig.suptitle('Numerical Quantum Mechanics: Matrix Eigenvalue & Basis Expansion Methods',
                 fontsize=14, fontweight='bold', y=0.99)

    #out = '/home/claude/ho_basis_results.png'
    #plt.savefig(out, dpi=150, bbox_inches='tight')
    #plt.close()
    #print(f"Figure saved to {out}")


if __name__ == "__main__":
    print_summary()
    plot_all()
    print("\nAll done.")
