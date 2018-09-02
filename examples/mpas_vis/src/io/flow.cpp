#include "flow.h"

#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkPolyLine.h>

flow::flow(int num_lines){

    flow::num_lines = num_lines;
    flow::fl_lengths.resize((num_lines));
}


flow::~flow(){

}

void flow::write_flow_to_vtk(std::string &filename){ // line, time-step, dimensions
    int n_steps = flow_lines.size()/(3*num_lines);



    vtkSmartPointer<vtkPoints> points =
            vtkSmartPointer<vtkPoints>::New();

    std::vector<double> p(3);
    int i = 0;
    // Iterate and print values of the list
    for (double n : flow_lines) {
        p[i%3] = n;
        i++;
        if (i%3==0){
            points->InsertNextPoint ( p[0],
                    p[1],
                    p[2]);
            //                std::cout<<"\n"<<p[0]<<" "<<p[1]<<" "<<p[2];

        }
    }
            std::cout<<" "<<i/3<<" "<<fl_lengths.size();"\n";


    // Create a cell array to store the lines in and add the lines to it
    vtkSmartPointer<vtkCellArray> cells =
            vtkSmartPointer<vtkCellArray>::New();
    int sum=0;
    for (int j=0;j<num_lines;j++){

        vtkSmartPointer<vtkPolyLine> polyLine =
                vtkSmartPointer<vtkPolyLine>::New();


//        polyLine->GetPointIds()->SetNumberOfIds(n_steps);
        polyLine->GetPointIds()->SetNumberOfIds(fl_lengths[j]);
//        for(unsigned int k = 0; k < n_steps; k++)
            for(unsigned int k = 0; k < fl_lengths[j]; k++)
        {

              polyLine->GetPointIds()->SetId(k,sum+k);
//            polyLine->GetPointIds()->SetId(k,j*n_steps+k);
        }

            if (fl_lengths[j]>1){
                cells->InsertNextCell(polyLine);
            }
                        sum += fl_lengths[j];

    }


    // Create a polydata to store everything in
    vtkSmartPointer<vtkPolyData> polydata =
            vtkSmartPointer<vtkPolyData>::New();

    // Add the points to the dataset
    polydata->SetPoints(points);

    // Add the lines to the dataset
    polydata->SetLines(cells);





    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
            vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(polydata);
#else
    writer->SetInputData(polydata);
#endif

    // Optional - set the mode. The default is binary.
    //writer->SetDataModeToBinary();
    //writer->SetDataModeToAscii();

    writer->Write();


}

//https://www.vtk.org/Wiki/VTK/Examples/Cxx/GeometricObjects/PolyLine
