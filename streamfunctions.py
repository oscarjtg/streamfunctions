import numpy as np

def streamfunction_direct_integration(u: np.ndarray, v: np.ndarray, xC: np.ndarray, xF: np.ndarray, yC: np.ndarray, yF: np.ndarray, psi_bot: np.ndarray, psi_left: np.ndarray) -> np.ndarray:
    """
    Directly integrates the definition of the stream function 
    using an explicit first order difference scheme.

    For stream function psi and velocity components u (x-direction) and v (y-direction)

    d(psi)/dy = u
    d(psi)/dx = -v

    where da/db is the partial derivative of a with respect to b.

    Parameters
    ----------
    u: np.ndarray 
        A two-dimensional NumPy array of shape `(m, n-1)` containing the 
        x-components of velocity in the X-Y plane on the vertical cell faces (xF, yC).

    v: np.ndarray 
        A two-dimensional NumPy array of shape `(m-1, n)` containing the 
        y-components of velocity in the X-Y plane on the horizontal cell faces (xC, yF).

    xC: np.ndarray
        A one-dimensional NumPy array of length `(m-1)` containing the x-coordinates 
        at the cell centres.

    xF: np.ndarray
        A one-dimensional NumPy array of length `(m)` containing the x-coordinates 
        at the cell faces.
    
    yC: np.ndarray
        A one-dimensional NumPy array of length `(n-1)` containing the y-coordinates 
        at the cell centres.

    yF: np.ndarray
        A one-dimensional NumPy array of length `(n)` containing the y-coordinates 
        at the cell faces.

    psi_bot: np.ndarray
        A one-dimensional NumPy array of length `(m)` containing the value of 
        the stream function along the bottom cells.

    psi_left: np.ndarray
        A one-dimensional NumPy array of length `(n)` containing the value of 
        the stream function at the cell corners along the left hand side of the domain.

    Returns
    -------
    np.ndarray 
        A two-dimensional NumPy array of shape `(m, n)` containing the value of 
        the computed stream function at the cell corners (xF, yF).
    """
    # Check that the inputs all have the correct dimensions.
    m = len(xF); n = len(yF)
    assert(u.shape == (m, n-1))
    assert(v.shape == (m-1, n))
    assert(len(xC) == m - 1)
    assert(len(yC) == n - 1)
    assert(len(psi_bot) == m)
    assert(len(psi_left) == n)
    assert(psi_bot[0] == psi_left[0])

    psi = np.zeros((m, n))
    psi[:, 0] = psi_bot
    for k in range(n - 1):
        psi[1:-1, k+1] = psi[1:-1, k] + u[1:-1, k] * (yF[k+1] - yF[k])

    psi[0, :] = psi_left
    for k in range(m - 1):
        psi[k+1, 1:-1] = psi[k, 1:-1] - v[k, 1:-1] * (xF[k+1] - xF[k])

    return psi

def vorticity(u: np.ndarray, v: np.ndarray, xC: np.ndarray, xF: np.ndarray, yC: np.ndarray, yF: np.ndarray) -> np.ndarray:
    """
    Calculates the vorticity of the 2D flow at cell corners (xF, yF).
    The vorticity is the curl of the velocity field.

    vort = dv/dx - du/dy

    Parameters
    ----------
    u: np.ndarray 
        A two-dimensional NumPy array of shape `(m, n-1)` containing the 
        x-components of velocity in the X-Y plane on the vertical cell faces (xF, yC).

    v: np.ndarray 
        A two-dimensional NumPy array of shape `(m-1, n)` containing the 
        y-components of velocity in the X-Y plane on the horizontal cell faces (xC, yF).

    xC: np.ndarray
        A one-dimensional NumPy array of length `(m-1)` containing the x-coordinates 
        at the cell centres.

    xF: np.ndarray
        A one-dimensional NumPy array of length `(m)` containing the x-coordinates 
        at the cell faces.
    
    yC: np.ndarray
        A one-dimensional NumPy array of length `(n-1)` containing the y-coordinates 
        at the cell centres.

    yF: np.ndarray
        A one-dimensional NumPy array of length `(n)` containing the y-coordinates 
        at the cell faces.

    Returns
    -------
    np.ndarray 
        A two-dimensional NumPy array of shape `(m, n)` containing the value of 
        the computed vorticity at the cell corners (xF, yF).
    """
    # Check that the inputs all have the correct dimensions.
    m = len(xF); n = len(yF)
    assert(u.shape == (m, n-1))
    assert(v.shape == (m-1, n))
    assert(len(xC) == m - 1)
    assert(len(yC) == n - 1)

    vort = np.zeros((m, n))
    dvdx = np.zeros((m, n))
    dudy = np.zeros((m, n))

    # This logic handles sides and corners elegantly on staggered grid.
    dvdx[1:-1, :] = (v[1:, :] - v[:-1, :]) / (xC[1:, np.newaxis] - xC[:-1, np.newaxis])
    dudy[:, 1:-1] = (u[:, 1:] - u[:, :-1]) / (yC[np.newaxis, 1:] - yC[np.newaxis, :-1])
    
    vort = dvdx - dudy

    return vort