// This is a script to get the TGeoManager from an edep-sim ROOT file OR GDML and extract the global position of a volume in a geometry. 
// To get the position,we need to recursively traverse the geometry tree until we encounter volume, taking account of the translations and rotations along the way.
// Author: Sam Fogarty, samuel.fogarty@colostate.edu

#include <TGeoManager.h>
#include <TFile.h>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <functional>

std::pair<bool, std::vector<std::vector<double>>> get_volume_position(const char* fileName, bool isROOT, const char* volName) {
    
    // Import the TGeoManager
    TGeoManager* geoManager = TGeoManager::Import(fileName);

    // Get the top node of the geometry
    TGeoNode* topVolume = geoManager->GetTopNode();

    // The global origins of the target volumes, to be determined below
    std::vector<std::vector<double>> globalOrigins;

    // Define a recursive function to traverse the geometry tree
    std::function<void(TGeoNode*, std::unique_ptr<TGeoHMatrix>)> traverse;

    traverse = [&traverse, &globalOrigins, &isROOT, &volName](TGeoNode* node, std::unique_ptr<TGeoHMatrix> globalMatrix) {
        if (globalMatrix == nullptr) {
            globalMatrix = std::make_unique<TGeoHMatrix>();  // identity matrix
        } else {
            // Combine the current node's transformation with the accumulated transformation
            TGeoHMatrix localMatrix = *globalMatrix;
            localMatrix.Multiply(node->GetMatrix());
            globalMatrix = std::make_unique<TGeoHMatrix>(localMatrix);
        }
        
        std::string suffix = "";
        if (isROOT) {
                suffix = "_PV";
        }
        
        // The volume names we're looking for. Note that right now this looks for one
        // volume at a time, but it does support finding multiple in a row.
        std::vector<std::string> targetVolumes = {volName + suffix};

        // Check if the node's volume is one of the volumes we're looking for
        for (const auto& volumeName : targetVolumes) {
            if (volumeName == node->GetVolume()->GetName()) { 
                Double_t localOrigin[3] = {0., 0., 0.};
                Double_t globalOrigin[3];
                // Convert local coordinates to global coordinates
                globalMatrix->LocalToMaster(localOrigin, globalOrigin);
                // Add the global origin to the list
                globalOrigins.push_back({globalOrigin[0], globalOrigin[1], globalOrigin[2]});
            }
        }

        // Recursively traverse the node's daughters
        for (int i = 0; i < node->GetNdaughters(); ++i) {
            TGeoNode* daughter = node->GetDaughter(i);
            traverse(daughter, std::make_unique<TGeoHMatrix>(*globalMatrix));
        }
    };

    // Start the traversal from the top volume
    traverse(topVolume, nullptr);
    
    return {true, globalOrigins};
}
