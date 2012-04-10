#include "common.hh"

/* ************************************************************************** *
 * NSSolver : Public methods                                                  *
 * ************************************************************************** */

NSSolver::NSSolver(u32 res, Scalar density, Scalar viscosity, Scalar pressure, Scalar diffusion) :
        mDelta(0), mTime(0), mRes(res), mDensity(density), mViscosity(viscosity),
                mPressure(pressure), mDiffusion(diffusion)
{
    mD = alloc();
    mOldD = alloc();
    mU = alloc();
    mOldU = alloc();
    mV = alloc();
    mOldV = alloc();
    mCurl = alloc();

    reset();
}

NSSolver::~NSSolver()
{
    free(mD);
    free(mOldD);
    free(mU);
    free(mOldU);
    free(mV);
    free(mOldV);
    free(mCurl);
}

void NSSolver::reset()
{
    fill(mD, mDensity);
    fill(mOldD, mDensity);
    fill(mU, 0);
    fill(mOldU, 0);
    fill(mV, 0);
    fill(mOldV, 0);
    fill(mCurl, 0);
}

void NSSolver::simulateStep(f64 delta)
{
    mTime += delta;
    mDelta = delta;

    velocityStep();
    densityStep();
}

/* ************************************************************************** *
 * NSSolver : Private methods                                                 *
 * ************************************************************************** */

NSSolver::Scalar** NSSolver::alloc()
{
    Scalar** array = new Scalar*[mRes + 2];
    for (u32 i = 0; i < mRes + 2; i++)
        array[i] = new Scalar[mRes + 2];

    return array;
}

void NSSolver::free(Scalar** array)
{
    for (u32 i = 0; i < mRes + 2; i++)
        delete array[i];

    delete array;
}

void NSSolver::fill(Scalar** array, Scalar value)
{
    for (u32 i = 0; i < mRes + 2; i++)
        std::fill(array[i], array[i] + mRes + 2, value);
}

void NSSolver::velocityStep()
{
    // add velocity that was input by mouse
    addSource(mU, mOldU);
    addSource(mV, mOldV);

    // add in vorticity confinement force
    vorticityConfinement();

    addSource(mU, mOldU);
    addSource(mV, mOldV);

    // add in buoyancy force
    buoyancy();
    addSource(mV, mOldV);

    // swapping arrays for economical mem use and calculating diffusion in velocity.
    std::swap(mU, mOldU);
    std::swap(mV, mOldV);

    diffuse(mU, mOldU, mViscosity);
    diffuse(mV, mOldV, mViscosity);

    // we create an incompressible field for more effective advection.
    project();
    std::swap(mU, mOldU);
    std::swap(mV, mOldV);

    // self advect velocities
    advect(1, mU, mOldU, mOldU, mOldV);
    advect(2, mV, mOldV, mOldU, mOldV);

    // make an incompressible field
    project();

    // clear all input velocities for next frame
    fill(mOldU, 0);
    fill(mOldV, 0);
}

void NSSolver::densityStep()
{
    // add density inputted by mouse
    addSource(mD, mOldD);
    std::swap(mD, mOldD);

    diffuse(mD, mOldD, mDiffusion);
    std::swap(mD, mOldD);

    advect(0, mD, mOldD, mU, mV);

    // clear input density array for next frame
    fill(mOldD, 0);
}

void NSSolver::addSource(Scalar** to, Scalar** src)
{
    for (u32 i = 0; i < mRes + 2; i++)
        for (u32 j = 0; j < mRes + 2; j++)
            to[i][j] += mDelta * src[i][j];
}

NSSolver::Scalar NSSolver::curl(u32 x, u32 y) const
{
    Scalar du_dx = (mU[x][y + 1] - mU[x][y - 1]); // * 0.5
    Scalar dv_dy = (mV[x + 1][x] - mV[x - 1][y]); // * 0.5

    return (du_dx - dv_dy) * 0.5;
}

void NSSolver::vorticityConfinement()
{
    // Calculate magnitude of curl(u,v) for each cell. (|w|)
    for (u32 i = 1; i <= mRes; i++)
        for (u32 j = 1; j <= mRes; j++)
            mCurl[i][j] = abs(curl(i, j));

#pragma omp parallel for
    for (u32 i = 2; i < mRes; i++)
        for (u32 j = 2; j < mRes; j++)
        {
            // Find derivative of the magnitude (n = del |w|)
            Scalar dw_dmU = (mCurl[i + 1][j] - mCurl[i - 1][j]) * 0.5f;
            Scalar dw_dmV = (mCurl[i][j + 1] - mCurl[i][j - 1]) * 0.5f;

            // Calculate vector length. (|n|)
            // Add small factor to prevent divide bmV zeros.
            Scalar length = sqrt(pow(dw_dmU, 2) + pow(dw_dmV, 2)) + 0.000001;

            // N = ( n/|n| )
            dw_dmU /= length;
            dw_dmV /= length;

            Scalar v = curl(i, j);

            // N mU w
            mOldU[i][j] = dw_dmV * -v;
            mOldV[i][j] = dw_dmU * v;
        }
}

void NSSolver::buoyancy()
{
    Scalar Tamb = 0;
    Scalar a = 0.000625f;
    Scalar b = 0.025f;

    // sum all temperatures
    for (u32 i = 1; i <= mRes; i++)
        for (u32 j = 1; j <= mRes; j++)
            Tamb += mD[i][j];

    // get average temperature
    Tamb /= pow(f64(mRes), 2);

    // for each cell compute buoyancy force
    for (u32 i = 1; i <= mRes; i++)
        for (u32 j = 1; j <= mRes; j++)
            mOldV[i][j] = a * mD[i][j] + -b * (mD[i][j] - Tamb);
}

void NSSolver::diffuse(Scalar** mU, Scalar** oldmU, Scalar diffusion)
{
    Scalar a = mDelta * diffusion * pow(mRes, 2);
    linearSolver(0, mU, oldmU, a, 1 + 4 * a);
}

void NSSolver::linearSolver(u32 b, Scalar** mU, Scalar** mU0, Scalar a, Scalar c)
{
    for (u32 k = 0; k < 20; k++)
    {
#pragma omp parallel for collapse(2)
        for (u32 i = 1; i <= mRes; i++)
            for (u32 j = 1; j <= mRes; j++)
                mU[i][j] = (a * (mU[i - 1][j] + mU[i + 1][j] + mU[i][j - 1] + mU[i][j + 1])
                        + mU0[i][j]) / c;

        setBoundaries(b, mU);
    }
}

void NSSolver::setBoundaries(u32 b, Scalar** mU)
{
    for (u32 i = 1; i <= mRes; i++)
    {
        Scalar v = (b == 1 ? -1 : 1);
        mU[0][i] = v * mU[1][i];
        mU[mRes + 1][i] = v * mU[mRes][i];

        v = (b == 2 ? -1 : 1);
        mU[i][0] = v * mU[i][1];
        mU[i][mRes + 1] = v * mU[i][mRes];
    }

    mU[0][0] = 0.5 * (mU[1][0] + mU[0][1]);
    mU[0][mRes + 1] = 0.5 * (mU[1][mRes + 1] + mU[0][mRes]);
    mU[mRes + 1][0] = 0.5 * (mU[mRes][0] + mU[mRes + 1][1]);
    mU[mRes + 1][mRes + 1] = 0.5 * (mU[mRes][mRes + 1] + mU[mRes + 1][mRes]);
}

void NSSolver::project()
{
    for (u32 i = 1; i <= mRes; i++)
        for (u32 j = 1; j <= mRes; j++)
        {
            mOldV[i][j] = (mU[i + 1][j] - mU[i - 1][j] + mV[i][j + 1] - mV[i][j - 1]) * -0.5 / mRes;
            mOldU[i][j] = 0;
        }

    setBoundaries(0, mOldV);
    setBoundaries(0, mOldU);

    linearSolver(0, mOldU, mOldV, 1, 4);

    for (u32 i = 1; i <= mRes; i++)
        for (u32 j = 1; j <= mRes; j++)
        {
            mU[i][j] -= 0.5f * mRes * (mOldU[i + 1][j] - mOldU[i - 1][j]);
            mV[i][j] -= 0.5f * mRes * (mOldU[i][j + 1] - mOldU[i][j - 1]);
        }

    setBoundaries(1, mU);
    setBoundaries(2, mV);
}

void NSSolver::advect(u32 b, Scalar** d, Scalar** d0, Scalar** du, Scalar** dv)
{
    u32 i0, j0, i1, j1;
    Scalar x, y, s0, t0, s1, t1, dt0;

    dt0 = mDelta * mRes;

    for (u32 i = 1; i <= mRes; i++)
        for (u32 j = 1; j <= mRes; j++)
        {
            // go backwards through velocity field
            x = i - dt0 * du[i][j];
            y = j - dt0 * dv[i][j];

            // interpolate results
            if (x > mRes + 0.5)
                x = mRes + 0.5f;
            if (x < 0.5)
                x = 0.5f;

            i0 = x;
            i1 = i0 + 1;

            if (y > mRes + 0.5)
                y = mRes + 0.5f;
            if (y < 0.5)
                y = 0.5f;

            j0 = y;
            j1 = j0 + 1;

            s1 = x - i0;
            s0 = 1 - s1;
            t1 = y - j0;
            t0 = 1 - t1;

            d[i][j] = s0 * (t0 * d0[i0][j0] + t1 * d0[i0][j1])
                    + s1 * (t0 * d0[i1][j0] + t1 * d0[i1][j1]);

        }

    setBoundaries(b, d);
}
