/*
 *  Zhiyuan Liu
 * Oct 22, 2018
 * This file is providing CLI for computing linking structure
 *
 */
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "ImageDistanceMap.h"
#include "RAWImageFile.h"
#include "M3DObjectFile.h"
#include "WorldSystem.h"
#include "AllImageIO.h"
#include <time.h>

M3DObject* readModel(const char * fileName)
{
        M3DObjectFile objectFile;
        SimilarityTransform3D * xform;
        M3DObject * obj;
        WorldSystem * world;

        if (fileName == NULL || strlen(fileName) == 0)
                return NULL;

        xform = new SimilarityTransform3D;
        world = new WorldSystem;
    int markedPrimitiveId = 0;
        obj = objectFile.read(fileName, markedPrimitiveId,
                xform, world, NULL);

        // cout << "Model file read successfully" << endl ;

        if (obj != NULL) {
        obj->setWorld(world, true);	// Sets world in object and object->that
        }
        else {
        delete xform;
                delete world;
                return NULL;
        }

        if (xform->wasRead()) {
#ifdef DEBUG
                cout << "Similarity transformation loaded" << endl;
#endif
        }
        else {			// The transformation in object is set to identity
                cout << "Similarity transformation not loaded" << endl;
        }
    obj->setTransformation(xform);
        obj->loadedObject()->setTransformation(xform);

        return obj;
}

int main(int argc, char** argv)
{
    if(argc < 3)
    {
        std::cout << "[Usage:] computeLinking [srepPath] [outputPath]" << std::endl;
        return -1;
    }

    char* modelName = argv[1];
    char* outputPath = argv[2];
    M3DObject* object = readModel(modelName);
    object->setOutputLinkingPath(outputPath);
    object->renderLinkingStructure();
}
