/******* COPYRIGHT ************************************************
*                                                                 *
*                             FlowVR                              *
*                     Application Library                         *
*                                                                 *
*-----------------------------------------------------------------*
 * COPYRIGHT (C) 2003-2011                by                       *
* INRIA                                                           *
* ALL RIGHTS RESERVED.	                                          *
*                                                                 *
* This source is covered by the GNU LGPL, please refer to the     *
* COPYING file for further information.                           *
*                                                                 *
*-----------------------------------------------------------------*
*                                                                 *
*  Original Contributors:                                         *
*    Jeremie Allard,                                              *
*    Thomas Arcila,                                               *
*    Jean-Denis Lesage.                                           *	
*    Clement Menier,                                              *
*    Bruno Raffin                                                 *
*                                                                 *
*******************************************************************
*                                                                 *
* File: get.cpp                                                   *
*                                                                 *
* Contacts:                                                       *
*  26/02/2004 Jeremie Allard <Jeremie.Allard@imag.fr>             *
*                                                                 *
******************************************************************/
#include "flowvr/module.h"
#include <iostream>
#include <unistd.h>

//Decaf
#include <decaf/data_model/simpleconstructdata.hpp>
#include <decaf/data_model/constructtype.h>
#include <decaf/data_model/arrayconstructdata.hpp>
#include <decaf/data_model/baseconstructdata.hpp>

using namespace decaf;
using namespace std;

int sleep_time=1;

int main(int argc, const char** argv)
{
  flowvr::InputPort pIn("text");
  std::vector<flowvr::Port*> ports;
  ports.push_back(&pIn);

  flowvr::ModuleAPI* flowvr = flowvr::initModule(ports);
  if (flowvr == NULL)
  {
    return 1;
  }

  int it=0;
  while (flowvr->wait())
  {
    // Get Message
    flowvr::Message m;
    flowvr->get(&pIn,m);

    //Creating the
    char* serialBuffer = (char*)(m.data.readAccess());
    int sizeSerialBuffer = m.data.getSize();
    std::shared_ptr<ConstructData> container = std::make_shared<ConstructData>();
    container->allocate_serial_buffer(sizeSerialBuffer);
    memcpy(container->getInSerialBuffer(), serialBuffer, sizeSerialBuffer);
    container->merge();

    //Extracting the fields
    std::shared_ptr<BaseConstructData> baseItData = container->getData("it");
    if(!baseItData)
    {
        std::cerr<<"ERROR : Could not find the field \"it\" in the datamodel. Abording."<<std::endl;
        exit(1);
    }
    std::shared_ptr<SimpleConstructData<int> > itData =
            dynamic_pointer_cast<SimpleConstructData<int> >(baseItData);

    std::shared_ptr<BaseConstructData> baseArrayData = container->getData("array");
    if(!baseArrayData)
    {
        std::cerr<<"ERROR : Could not find the field \"array\" in the datamodel. Abording."<<std::endl;
        exit(1);
    }
    std::shared_ptr<ArrayConstructData<int> > arrayData =
            dynamic_pointer_cast<ArrayConstructData<int> >(baseArrayData);
    int* array = arrayData->getArray();


    std::cout<<"Received : "<<itData->getData()<<std::endl;
    std::cout<<"[";
    for(unsigned int i = 0; i < 10; i++)
        std::cout<<array[i]<<" ";
    std::cout<<"]"<<std::endl;

    sleep(sleep_time);
    ++it;
  }

  flowvr->close();
  return 0;
}