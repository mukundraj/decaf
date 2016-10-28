//---------------------------------------------------------------------------
//
// 2-node producer-consumer coupling example
//
// prod (4 procs) - con (2 procs)
//
// entire workflow takes 8 procs (2 dataflow procs between prod and con)
// this file contains the consumer (2 procs)
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// mdreher@anl.gov
//
//--------------------------------------------------------------------------

#include <decaf/decaf.hpp>
#include <decaf/data_model/simplefield.hpp>
#include <decaf/data_model/arrayfield.hpp>
#include <decaf/data_model/blockfield.hpp>
#include <decaf/data_model/array3dconstructdata.hpp>
#include <decaf/data_model/boost_macros.h>

#include <decaf/data_model/morton.h>

#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <map>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <set>

#include "Damaris.h"

#include <decaf/workflow.hpp>

using namespace decaf;
using namespace std;
using namespace boost;


enum targetType {
    ABS = 0,
    REL
};

#define MAX_SIZE_REQUEST 2048
typedef struct
{
    int type;                   //Target absolute in space or relative to a reference point
    float target[3];            //Coordonates of the absolute target
    char targetRequest[2048];   //Request of the
}Target;

float length(float* vec){
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

void normalize(float* vec){
    float l = length(vec);
    vec[0] = vec[0] / l;
    vec[1] = vec[1] / l;
    vec[2] = vec[2] / l;
}

float distance(float* p1, float* p2)
{
    float vec[3];
    vec[0] = p2[0] - p1[0];
    vec[1] = p2[1] - p1[1];
    vec[2] = p2[2] - p1[2];
    return length(vec);
}

std::set<int> filterIds;
std::vector<int> arrayIds; //HARD CODED for SimpleWater example

std::vector<Target> targets;

float distanceValidTarget = 0.5;
float maxTotalForces = 3.f;
std::string model;

void loadTargets()
{
    if(model.compare(std::string("SimplePeptideWater")) == 0)
    {
        Target target;
        target.target[0] = 4.0;
        target.target[1] = 6.0;
        target.target[2] = 20.0;
        targets.push_back(target);

        target.target[0] = 23.0;
        target.target[1] = 6.0;
        target.target[2] = 20.0;
        targets.push_back(target);

        target.target[0] = 23.0;
        target.target[1] = 38.0;
        target.target[2] = 20.0;
        targets.push_back(target);

        target.target[0] = 4.0;
        target.target[1] = 38.0;
        target.target[2] = 20.0;
        targets.push_back(target);

        filterIds = {  109, 110, 111, 112, 113, 114,
                       115, 116, 117, 118, 119, 120
                    }; //HARD CODED for SimpleWater example

        arrayIds = {  109, 110, 111, 112, 113, 114,
                      115, 116, 117, 118, 119, 120
                   }; //HARD CODED for SimpleWater example
    }
    else if(model.compare(std::string("fepa")) == 0)
    {
        Target target;
        target.target[0] = 50.649998;
        target.target[1] = 40.020000;
        target.target[2] = 74.940002;
        targets.push_back(target);

        target.target[0] = 57.994247;
        target.target[1] = 42.744064;
        target.target[2] = 75.205559;
        targets.push_back(target);

        target.target[0] = 58.028599;
        target.target[1] = 39.480324;
        target.target[2] = 62.716755;
        targets.push_back(target);

        target.target[0] = 58.175446;
        target.target[1] = 36.721069;
        target.target[2] = 59.135941;
        targets.push_back(target);

        target.target[0] = 60.568310;
        target.target[1] = 35.987762;
        target.target[2] = 56.373985;
        targets.push_back(target);

        target.target[0] = 57.443069;
        target.target[1] = 41.200779;
        target.target[2] = 52.448627;
        targets.push_back(target);

        target.target[0] = 60.272179;
        target.target[1] = 41.596397;
        target.target[2] = 41.934307;
        targets.push_back(target);

        target.target[0] = 58.013557;
        target.target[1] = 49.347263;
        target.target[2] = 14.191130;
        targets.push_back(target);

        //Ids for ENT and FE residues
        for(int i = 69901; i <= 69952; i++)
        {
            filterIds.insert(i);
            arrayIds.push_back(i);
        }
    }
}


// consumer
void target(Decaf* decaf)
{
    vector< pConstructData > in_data;
    fprintf(stderr, "Launching target\n");
    fflush(stderr);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int iteration = 0;

    loadTargets();

    unsigned int currentTarget = 0;

    while (decaf->get(in_data))
    {
        ArrayFieldf posField = in_data[0]->getFieldData<ArrayFieldf>("pos");
        if(!posField)
        {
            fprintf(stderr, "ERROR : the field \'pos\' is not in the data model\n");
            continue;
        }

        float* pos = posField.getArray();

        ArrayFieldu idsField = in_data[0]->getFieldData<ArrayFieldu>("ids");
        if(!idsField)
        {
            fprintf(stderr, "ERROR : the field \'ids\' is not in the data model\n");
            continue;
        }

        std::shared_ptr<Array3DConstructData<float> > grid =
                in_data[0]->getTypedData<Array3DConstructData<float> >("grid");

        unsigned int* ids = idsField.getArray();

        unsigned int nbAtoms = posField->getNbItems();

        float avg[3];
        bzero(avg, 3 * sizeof(float));
        int count = 0;

        // Filtering the data to select only the steered complex
        for(unsigned int i = 0; i < nbAtoms; i++)
        {
            if(filterIds.count(ids[i]) != 0)
            {
                avg[0] = avg[0] + pos[3*i];
                avg[1] = avg[1] + pos[3*i+1];
                avg[2] = avg[2] + pos[3*i+2];
                count++;
            }
        }

        // Computing the average position of the steered system
        avg[0] = avg[0] / (float)filterIds.size();
        avg[1] = avg[1] / (float)filterIds.size();
        avg[2] = avg[2] / (float)filterIds.size();

        fprintf(stderr, "Average position : %f %f %f\n", avg[0], avg[1], avg[2]);

        // Checking the position between the system and the target
        float* targetPos = targets[currentTarget].target;
        float dist = distance(avg, targetPos);
        if(dist < distanceValidTarget)
        {
            //Changing the active target
            currentTarget++;
            fprintf(stderr, "Changing target\n");
            if(currentTarget < targets.size())
                targetPos = targets[currentTarget].target;
            else
            {
                fprintf(stderr, "Steering terminated. Closing the app\n");
                break;
            }
        }
        fprintf(stderr, "[%i/%u] Target position: %f %f %f\n", currentTarget, targets.size(), targetPos[0], targetPos[1], targetPos[2]);
        fprintf(stderr, "[%i/%u] Distance to target: %f\n", currentTarget, targets.size(), dist);

        float force[3];

        if(grid)
        {
            fprintf(stderr, "Grid available, we are using path finding.\n");
            // Computing the force direction
            force[0] = targetPos[0] - avg[0];
            force[1] = targetPos[1] - avg[1];
            force[2] = targetPos[2] - avg[2];

            normalize(force);

            // Computing the force intensity
            float maxForcesPerAtom = maxTotalForces / (float)filterIds.size();
            force[0] = force[0] * maxForcesPerAtom;
            force[1] = force[1] * maxForcesPerAtom;
            force[2] = force[2] * maxForcesPerAtom;

            //Block<3> block = grid->getBlock();

            // Getting the grid info
            BlockField blockField  = in_data[0]->getFieldData<BlockField>("domain_block");
            if(!blockField)
            {
                fprintf(stderr, "ERROR: the field \'domain_block\' is not in the data model\n");
                continue;
            }
            Block<3>* domainBlock = blockField.getBlock();

            // We use the domain because we gather the full grid on 1 node

            // Using the local coords because the 3D array is a local
            // representation, not own representation
            unsigned int avgIndex[3];
            if(!domainBlock->getLocalPositionIndex(avg, avgIndex))
            {
                fprintf(stderr," ERROR: Unable to get the correct index of the tracked position\n");
                continue;
            }
            unsigned int targetIndex[3];
            if(!domainBlock->getLocalPositionIndex(targetPos,targetIndex))
            {
                fprintf(stderr," ERROR: Unable to get the correct index of the target position\n");
                continue;
            }

            fprintf(stderr,"Track cell: [%u %u %u]\n", avgIndex[0], avgIndex[1], avgIndex[2]);
            fprintf(stderr,"Target cell: [%u %u %u]\n", targetIndex[0], targetIndex[1], targetIndex[2]);

            // Now we can do the path finding


        }
        else
        {
            fprintf(stderr, "Computing direct force\n");
            // Computing the force direction
            force[0] = targetPos[0] - avg[0];
            force[1] = targetPos[1] - avg[1];
            force[2] = targetPos[2] - avg[2];

            normalize(force);

            // Computing the force intensity
            float maxForcesPerAtom = maxTotalForces / (float)filterIds.size();
            force[0] = force[0] * maxForcesPerAtom;
            force[1] = force[1] * maxForcesPerAtom;
            force[2] = force[2] * maxForcesPerAtom;
        }


        fprintf(stderr, "Force emitted : %f %f %f\n", force[0], force[1], force[2]);

        // Generating the data model
        pConstructData container;
        ArrayFieldf field(force, 3, 3, false);
        container->appendData("force", field, DECAF_NOFLAG, DECAF_PRIVATE,
                              DECAF_SPLIT_KEEP_VALUE, DECAF_MERGE_FIRST_VALUE);
        if(iteration == 0)
        {
            ArrayFieldi fieldIds(&arrayIds[0], arrayIds.size(), 1, false);
            container->appendData("ids", fieldIds, DECAF_NOFLAG, DECAF_PRIVATE,
                                  DECAF_SPLIT_DEFAULT, DECAF_MERGE_APPEND_VALUES);
        }

        decaf->put(container);

        //Clearing the first time as the data model will change at the next iteration
        if(iteration == 0)
            decaf->clearBuffers(DECAF_NODE);

        iteration++;
    }

    // terminate the task (mandatory) by sending a quit message to the rest of the workflow
    fprintf(stderr, "Target terminating\n");
    decaf->terminate();
}

// every user application needs to implement the following run function with this signature
// run(Workflow&) in the global namespace
void run(Workflow& workflow)                             // workflow
{
    MPI_Init(NULL, NULL);

    char processorName[MPI_MAX_PROCESSOR_NAME];
    int size_world, rank, nameLen;

    MPI_Comm_size(MPI_COMM_WORLD, &size_world);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processorName,&nameLen);

    srand(time(NULL) + rank * size_world + nameLen);

    fprintf(stderr, "target rank %i\n", rank);

    // create decaf
    Decaf* decaf = new Decaf(MPI_COMM_WORLD, workflow);


    // start the task
    target(decaf);

    // cleanup
    delete decaf;
    fprintf(stderr,"Decaf deleted. Waiting on finalize\n");
    MPI_Finalize();
}

// test driver for debugging purposes
// normal entry point is run(), called by python
int main(int argc,
         char** argv)
{
    fprintf(stderr, "Hello treatment\n");
    // define the workflow
    Workflow workflow;
    //make_wflow(workflow);
    Workflow::make_wflow_from_json(workflow, "wflow_gromacs.json");

    if(argc != 4)
    {
        fprintf(stderr, "Usage : targetmanager profile maxforce disttotarget\n");
        return 0;
    }

    model = string(argv[1]);
    maxTotalForces = atof(argv[2]);
    distanceValidTarget = atof(argv[3]);

    // run decaf
    run(workflow);

    return 0;
}
