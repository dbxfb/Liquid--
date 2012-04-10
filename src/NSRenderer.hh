#ifndef NSRENDERER_HH_
#define NSRENDERER_HH_

class NSRenderer
{
public:
    NSRenderer(int* argc, char** argv, u32 size, NSSolver& solver);
    ~NSRenderer();

    void run();
    void post()
    {
        mHasData = true;
    }

private:
    static void doMenu(int menu);
    static void doRefresh();
    static void doRender();
    static void doMouse(int button, int state, int x, int y);
    static void doMouseMove(int x, int y);

    void createMenu();
    void createVertices();
    void createColors();

    void computeStep();
    void render();

    void handleMenu(int menu);
    void mouseClick(u32 button, u32 state, u32 x, u32 y);
    void mouseMove(u32 x, u32 y);

    u32 mSize;
    NSSolver& mSolver;
    bool mHasData;

    f32* mVertices;
    u8* mColors;

    u32 mMouseButton;
    bool mMouseDown;
    s32 mPrevX, mPrevY;

    bool mPixelBlending;
    bool mVelocityVectors;

    static NSRenderer* mInstance;
};

#endif /* NSRENDERER_HH_ */
