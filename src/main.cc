#include "common.hh"

using namespace std;

static NSSolver* Solver = NULL;
static NSRenderer* Renderer = NULL;

void* computationsThread(void*)
{
    std::queue<f64> timings;
    f64 total = 0;
    u32 i = 1;

    Timer t, bench;
    while (true)
    {
        Solver->simulateStep(t.delta());

        f64 time = bench.delta();
        total += time;
        timings.push(time);
        if(timings.size() > 13)
        {
            total -= timings.front();
            timings.pop();
            i++;
        }

        if(i % 13 == 0)
            cout << (13 / total) << " fps" << endl;

        Renderer->post();
    }

    return NULL;
}

int main(int argc, char** argv)
{
    Solver = new NSSolver(250, 0, 0, 1, 0);
    Renderer = new NSRenderer(&argc, argv, 750, *Solver);

    pthread_t computations;
    pthread_create(&computations, NULL, &computationsThread, NULL);

    Renderer->run();

    return 0;
}


