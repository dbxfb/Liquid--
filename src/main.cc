#include "common.hh"

using namespace std;

static NSSolver* Solver = NULL;
static NSRenderer* Renderer = NULL;

void* computationsThread(void*)
{
    Timer t, bench;
    while (true)
    {
        Solver->simulateStep(t.delta());
        cout << (bench.delta() * 1000) << endl;
        Renderer->post();
    }

    return NULL;
}

int main(int argc, char** argv)
{
    Solver = new NSSolver(256, 1, 0, 1, 0);
    Renderer = new NSRenderer(&argc, argv, 1024, *Solver);

    pthread_t computations;
    pthread_create(&computations, NULL, &computationsThread, NULL);

    Renderer->run();

    return 0;
}


