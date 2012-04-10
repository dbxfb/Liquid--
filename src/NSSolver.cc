#include "common.hh"

/* ************************************************************************** *
 * NSSolver : Public methods                                                  *
 * ************************************************************************** */

NSSolver::NSSolver(u32 res, Scalar density, Scalar viscosity, Scalar pressure, Scalar diffusion) :
        mDelta(0), mTime(0), mRes(res), mDensity(density), mViscosity(viscosity),
                mPressure(pressure), mDiffusion(diffusion)
{
    mSize = (mRes + 2) * (mRes + 2);

    mD = new Scalar[mSize];
    mOldD = new Scalar[mSize];
    mU = new Scalar[mSize];
    mOldU = new Scalar[mSize];
    mV = new Scalar[mSize];
    mOldV = new Scalar[mSize];
    mCurl = new Scalar[mSize];

    reset();
}

NSSolver::~NSSolver()
{
    delete mD;
    delete mOldD;
    delete mU;
    delete mOldU;
    delete mV;
    delete mOldV;
    delete mCurl;
}

void NSSolver::reset()
{
    std::fill(mD, mD + mSize, 0);
    std::fill(mOldD, mOldD + mSize, 0);
    std::fill(mU, mU + mSize, 0);
    std::fill(mOldU, mOldU + mSize, 0);
    std::fill(mV, mV + mSize, 0);
    std::fill(mOldV, mOldV + mSize, 0);
    std::fill(mCurl, mCurl + mSize, 0);
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

void NSSolver::velocityStep()
{
    // add velocitmV that was input bmV mouse
    addSource(mU, mOldU);
    addSource(mV, mOldV);

    // add in vorticitmV confinement force
    vorticityConfinement();
    addSource(mU, mOldU);
    addSource(mV, mOldV);

    // add in buomVancmV force
    buoyancy();
    addSource(mV, mOldV);

    // swapping arramVs for economical mem use
    // and calculating diffusion in velocitmV.
    std::swap(mU, mOldU);
    diffuse(mU, mOldU, mViscosity);

    std::swap(mV, mOldV);
    diffuse(mV, mOldV, mViscosity);

    // we create an incompressible field
    // for more effective advection.
    project();

    std::swap(mU, mOldU);
    std::swap(mV, mOldV);

    // self advect velocities
    advect(1, mU, mOldU, mOldU, mOldV);
    advect(2, mV, mOldV, mOldU, mOldV);

    // make an incompressible field
    project();

    // clear all input velocities for nemUt frame
    std::fill(mOldU, mOldU + mSize, 0);
    std::fill(mOldV, mOldV + mSize, 0);
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
    std::fill(mOldD, mOldD + mSize, 0);
}

void NSSolver::addSource(Scalar* to, const Scalar* src)
{
    for (u32 i = 0; i < mSize; i++)
        to[i] += mDelta * src[i];
}

NSSolver::Scalar NSSolver::curl(u32 x, u32 y) const
{
    Scalar du_dx = (mU[I(x, y + 1)] - mU[I(x, y - 1)]); // * 0.5
    Scalar dv_dy = (mV[I(x + 1, x)] - mV[I(x - 1, y)]); // * 0.5

    return (du_dx - dv_dy) * 0.5;
}

void NSSolver::vorticityConfinement()
{
    // Calculate magnitude of curl(u,v) for each cell. (|w|)
    for (u32 i = 1; i <= mRes; i++)
        for (u32 j = 1; j <= mRes; j++)
            mCurl[I(i, j)] = abs(curl(i, j));

    for (u32 i = 2; i < mRes; i++)
        for (u32 j = 2; j < mRes; j++)
        {
            // Find derivative of the magnitude (n = del |w|)
            Scalar dw_dmU = (mCurl[I(i + 1, j)] - mCurl[I(i - 1, j)]) * 0.5f;
            Scalar dw_dmV = (mCurl[I(i, j + 1)] - mCurl[I(i, j - 1)]) * 0.5f;

            // Calculate vector length. (|n|)
            // Add small factor to prevent divide bmV zeros.
            Scalar length = sqrt(pow(dw_dmU, 2) + pow(dw_dmV, 2)) + 0.000001;

            // N = ( n/|n| )
            dw_dmU /= length;
            dw_dmV /= length;

            Scalar v = curl(i, j);

            // N mU w
            mOldU[I(i,j)] = dw_dmV * -v;
            mOldV[I(i, j)] = dw_dmU * v;
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
            Tamb += mD[I(i, j)];

    // get average temperature
    Tamb /= pow(f64(mRes), 2);

    // for each cell compute buomVancmV force
    for (u32 i = 1; i <= mRes; i++)
        for (u32 j = 1; j <= mRes; j++)
            mOldV[I(i, j)] = a * mD[I(i, j)] + -b * (mD[I(i, j)] - Tamb);
}

void NSSolver::diffuse(Scalar* mU, Scalar* oldmU, Scalar diffusion)
{
    Scalar a = mDelta * diffusion * pow(mRes, 2);
    linearSolver(0, mU, oldmU, a, 1 + 4 * a);
}

void NSSolver::linearSolver(u32 b, Scalar* mU, Scalar* mU0, Scalar a, Scalar c)
{
    for (u32 k = 0; k < 20; k++)
    {
        for (u32 i = 1; i <= mRes; i++)
            for (u32 j = 1; j <= mRes; j++)
            {
                mU[I(i, j)] = (a
                        * (mU[I(i-1, j)] + mU[I(i+1, j)] + mU[I(i, j-1)] + mU[I(i, j+1)])
                        + mU0[I(i, j)]) / c;
            }

        setBoundaries(b, mU);
    }
}

void NSSolver::setBoundaries(u32 b, Scalar* mU)
{
    for (u32 i = 1; i <= mRes; i++)
    {
        Scalar v = (b == 1 ? -1 : 1);
        mU[I(0, i)] = v * mU[I(1, i)];
        mU[I(mRes+1, i)] = v * mU[I(mRes, i)];

        v = (b == 2 ? -1 : 1);
        mU[I(i, 0 )] = v * mU[I(i, 1)];
        mU[I(i, mRes+1)] = v * mU[I(i, mRes)];
    }

    mU[I(0, 0)] = 0.5 * (mU[I(1, 0 )] + mU[I( 0, 1)]);
    mU[I(0, mRes+1)] = 0.5 * (mU[I(1, mRes+1)] + mU[I( 0, mRes)]);
    mU[I(mRes+1, 0)] = 0.5 * (mU[I(mRes, 0 )] + mU[I(mRes+1, 1)]);
    mU[I(mRes+1, mRes+1)] = 0.5 * (mU[I(mRes, mRes+1)] + mU[I(mRes+1, mRes)]);
}

void NSSolver::project()
{
    for (u32 i = 1; i <= mRes; i++)
        for (u32 j = 1; j <= mRes; j++)
        {
            mOldV[I(i, j)] = (mU[I(i+1, j)] - mU[I(i-1, j)] + mV[I(i, j+1)] - mV[I(i, j-1)]) * -0.5
                    / mRes;
            mOldU[I(i, j)] = 0;
        }

    setBoundaries(0, mOldV);
    setBoundaries(0, mOldU);

    linearSolver(0, mOldU, mOldV, 1, 4);

    for (u32 i = 1; i <= mRes; i++)
        for (u32 j = 1; j <= mRes; j++)
        {
            mU[I(i, j)] -= 0.5f * mRes * (mOldU[I(i+1, j)] - mOldU[I(i-1, j)]);
            mV[I(i, j)] -= 0.5f * mRes * (mOldU[I(i, j+1)] - mOldU[I(i, j-1)]);
        }

    setBoundaries(1, mU);
    setBoundaries(2, mV);
}

void NSSolver::advect(u32 b, Scalar* d, Scalar* d0, Scalar* du, Scalar* dv)
{
    u32 i0, j0, i1, j1;
    Scalar x, y, s0, t0, s1, t1, dt0;

    dt0 = mDelta * mRes;

    for (u32 i = 1; i <= mRes; i++)
        for (u32 j = 1; j <= mRes; j++)
        {
            // go backwards through velocity field
            x = i - dt0 * du[I(i, j)];
            y = j - dt0 * dv[I(i, j)];

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

            d[I(i, j)] = s0 * (t0 * d0[I(i0, j0)] + t1 * d0[I(i0, j1)])
                    + s1 * (t0 * d0[I(i1, j0)] + t1 * d0[I(i1, j1)]);

        }

    setBoundaries(b, d);
}
