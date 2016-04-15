#include <OpenMesh/Core/IO/MeshIO.hh>
#include <kt84/util.hh>
#include <patchgen/generate_topology.hh>
#include "decl.hh"
#include "Patch.hh"
using namespace std;
using namespace Eigen;
using namespace kt84;

int main_glut(int argc, char* argv[]);

int main(int argc, char* argv[]) {
    int mode = util::simple_prompt<int>("mode (0=manually specify input, 1=random input, 2=GLUT app)");
    if (mode < 0 || 2 < mode) {
        cout << "unknown mode!\n";
        return 0;
    }
    
    if (mode == 1 || mode == 2) {
        unsigned long seed = util::simple_prompt<unsigned long>("random seed (0=current time)");
        if (seed == 0) {
            seed = static_cast<unsigned long>(std::time(nullptr));
            cout << "seed=" << seed << endl;
        }
        util::srand(seed);
    }
    
    if (mode == 2)
        return main_glut(argc, argv);
    
    while (true) {
        if (mode == 0) {
            VectorXi l;
            cout << "input number of edge subdivisions per side (finish by typing 0): ";
            while (true) {
                int i;
                cin >> i;
                if (i <= 0) break;
                l.conservativeResize(l.size() + 1);
                l[l.size() - 1] = i;
            }
            
            int num_sides = l.size();
            
            if (num_sides < 2 || 6 < num_sides) {
                cout << "num_sides=" << num_sides << " is unsupported.\n";
                continue;
            }
            if (l.sum() % 2 != 0) {
                cout << "the sum of number of edge subdivisions should be even.\n";
                continue;
            }
            if (l.sum() < 4) {
                cout << "input numbers are too small.\n";
                continue;
            }
            
            patchgen::PatchParam param;
            demo::Patch patch;
            patchgen::generate_topology(l, param, patch);
            demo::determine_geometry(patch, param.l);
            demo::save_mesh(patch, demo::get_fname(param));
        
        } else {
            int num_sides = util::simple_prompt<int>("num_sides");
            if (num_sides < 2 || 6 < num_sides) {
                cout << "num_sides=" << num_sides << " is unsupported.\n";
                continue;
            }
            int max_subdiv = util::simple_prompt<int>("max_subdiv");
            int num_tests  = util::simple_prompt<int>("num_tests");
            for (int k = 0; k < num_tests; ++k) {
                cout << "k=" << k << ",  ";
                VectorXi l = VectorXi::Zero(num_sides);
                for (int i = 0; i < num_sides; ++i)
                    l[i] = util::random_int(1, max_subdiv);
                if (l.sum() % 2 != 0) ++l[0];
                if (l.sum() < 4) l[0] += 2;
                
                patchgen::PatchParam param;
                demo::Patch patch;
                patchgen::generate_topology(l, param, patch);
                demo::determine_geometry(patch, param.l);
                demo::save_mesh(patch, demo::get_fname(param));
            }
        }
    }
    return 0;
}
