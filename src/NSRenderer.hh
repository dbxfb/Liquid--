#ifndef NSRENDERER_HH_
#define NSRENDERER_HH_

class NSRenderer
{
public:
    NSRenderer(int* argc, char** argv, u32 size, NSSolver& solver);
    ~NSRenderer();

    void post();
    void run();

private:
    static void doRender();
    static void doMouse(int button, int state, int x, int y);
    static void doMouseMove(int x, int y);

    void createVertices();
    void createColors();

    void render();

    void mouseClick(u32 button, u32 state, u32 x, u32 y);
    void mouseMove(u32 x, u32 y);

    void makeVec(f32 x, f32 y, f32 dx, f32 dy);

    u32 mSize;
    NSSolver& mSolver;

    f32* mVertices;
    u8* mColors;

    u32 mMouseButton;
    bool mMouseDown;
    s32 mPrevX, mPrevY;

    static NSRenderer* mInstance;
};


#endif /* NSRENDERER_HH_ */
