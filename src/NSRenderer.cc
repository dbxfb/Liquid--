#include "common.hh"

#define GL_GLEXT_PROTOTYPES 1

#include <GL/glut.h>

/* ************************************************************************** *
 * NSRenderer : Public methods                                                *
 * ************************************************************************** */

NSRenderer::NSRenderer(int* argc, char** argv, u32 size, NSSolver& solver) :
        mSize(size), mSolver(solver)
{
    mInstance = this;

    glutInit(argc, argv);
    glutInitWindowSize(mSize, mSize);
    glutInitWindowPosition(50, 50);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutCreateWindow("Smoke !");
    glClearColor(.3, .3, .5, 0);
    gluOrtho2D(0, mSize, 0, mSize);

    createVertices();
    createColors();
}

NSRenderer::~NSRenderer()
{
    delete mVertices;
    delete mColors;
}

void NSRenderer::run()
{
    glutDisplayFunc(&NSRenderer::doRender);
    glutMouseFunc(&NSRenderer::doMouse);
    glutMotionFunc(&NSRenderer::doMouseMove);
    glutIdleFunc(&NSRenderer::doRefresh);

    // Has to be run as the main program thread !
    glutMainLoop();
}

/* ************************************************************************** *
 * NSRenderer : Private methods                                               *
 * ************************************************************************** */
void NSRenderer::createVertices()
{
    u32 solverSize = mSolver.getRes();
    u32 line = solverSize * 8;

    mVertices = new f32[solverSize * line];

    f32 dg = f32(mSize) / f32(solverSize);
    f32 dg_2 = dg / 2.0;

    for (u32 i = 0; i < solverSize; i++)
    {
        u32 dx = (i + 0.5) * dg;

        for (u32 j = 0; j < solverSize; j++)
        {
            u32 dy = (j + 0.5) * dg;

            u32 x = dx - dg_2;
            u32 y = dy - dg_2;

            u32 idx = i * 8 + j * line;
            mVertices[idx + 0] = x;
            mVertices[idx + 1] = y;
            mVertices[idx + 2] = x + dg;
            mVertices[idx + 3] = y;
            mVertices[idx + 4] = x + dg;
            mVertices[idx + 5] = y + dg;
            mVertices[idx + 6] = x;
            mVertices[idx + 7] = y + dg;
        }
    }

}

void NSRenderer::createColors()
{
    u32 solverSize = mSolver.getRes();
    u32 line = solverSize * 4 * 3;

    mColors = new u8[solverSize * line];
}

void NSRenderer::render()
{
    u32 solverSize = mSolver.getRes();
    u32 line = solverSize * 4 * 3;

#if 0
    for (u32 i = 0; i < solverSize; i++)
    for (u32 j = 0; j < solverSize; j++)
    {
        u32 idx = i * 4 * 3 + j * line;
        f32 d = std::max(0.0, 1.0 - mSolver.getDensity(i + 1, j + 1));
        u32 c = d * 255;
        std::fill(mColors + idx, mColors + idx + 4 * 3, c);
    }
#else
    // Hardcore pixel blending from hell :)
    for (u32 j = 0; j < solverSize; j++)
    {
        u8* ptr = mColors + (j * line);
        for (u32 i = 0; i < solverSize; i++)
        {
            u8 d;
            {
                d = 255.0 * std::max(0.0, 1.0 - mSolver.getDensity(i, j));
                std::fill(ptr, ptr + 3, d);
                ptr += 3;
            }
            {
                d = 255.0 * std::max(0.0, 1.0 - mSolver.getDensity(i + 1, j));
                std::fill(ptr, ptr + 3, d);
                ptr += 3;
            }
            {
                d = 255.0 * std::max(0.0, 1.0 - mSolver.getDensity(i + 1, j + 1));
                std::fill(ptr, ptr + 3, d);
                ptr += 3;
            }
            {
                d = 255.0 * std::max(0.0, 1.0 - mSolver.getDensity(i, j + 1));
                std::fill(ptr, ptr + 3, d);
                ptr += 3;
            }
        }
    }
#endif

    glClear(GL_COLOR_BUFFER_BIT);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    {
        glVertexPointer(2, GL_FLOAT, 0, mVertices);
        glColorPointer(3, GL_UNSIGNED_BYTE, 0, mColors);
        glDrawArrays(GL_QUADS, 0, solverSize * solverSize * 4);
    }
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);

    glutSwapBuffers();
}

void NSRenderer::makeVec(f32 x, f32 y, f32 dx, f32 dy)
{
    glBegin(GL_LINES);
    glVertex2f(x, y);
    glVertex2f(x + dx, y + dy);
    glEnd();
}

void NSRenderer::mouseClick(u32 button, u32 state, u32 x, u32 y)
{
    mMouseButton = button;
    mMouseDown = (state == GLUT_DOWN);
}

void NSRenderer::mouseMove(u32 x, u32 y)
{
    if (!mMouseDown)
        return;

    u32 solverSize = mSolver.getRes();

    // get index for fluid cell under mouse position
    u32 i = (f64(x) / 1024.0) * solverSize + 1;
    u32 j = (f64(y) / 1024.0) * solverSize + 1;

    // set boundries
    if (i > solverSize)
        i = solverSize;
    if (i < 1)
        i = 1;
    if (j > solverSize)
        j = solverSize;
    if (j < 1)
        j = 1;

    j = solverSize - j;

    if (mMouseButton == GLUT_LEFT_BUTTON)
        mSolver.setDensity(i, j, 100);

    if (mMouseButton == GLUT_RIGHT_BUTTON)
    {
        s32 U = (s32(i) - mPrevX) * 5;
        s32 V = (s32(j) - mPrevY) * 5;
        mSolver.setU(i, j, U);
        mSolver.setV(i, j, V);
    }

    mPrevX = i;
    mPrevY = j;
}

/* ************************************************************************** *
 * NSRenderer : Static methods                                                *
 * ************************************************************************** */

NSRenderer* NSRenderer::mInstance = NULL;

void NSRenderer::doRefresh()
{
    if (mInstance->mHasData)
    {
        glutPostRedisplay();
        mInstance->mHasData = false;
    }
    else
        usleep(10000);
}

void NSRenderer::doRender()
{
    mInstance->render();
}

void NSRenderer::doMouse(int button, int state, int x, int y)
{
    mInstance->mouseClick(button, state, x, y);
}

void NSRenderer::doMouseMove(int x, int y)
{
    mInstance->mouseMove(x, y);
}

