import fire
import numpy as np
import ROOT
import csv

def get_volume_position(input_file, isROOT, volName):
    """
    Function that calls a ROOT cpp macro that calculates the module offsets from
    an edep-sim ROOT file OR geometry GDML.
    
    Will return an empty list if the active LAr volumes are not found.
    Note that the cpp script has a handful of hard-coded volume names in order to find the module
    offsets. So if the script isn't finding the offsets, you should check the cpp script and the
    GDML and make sure it's looking for the right volumes.
    
    Args:
        input_file (str): path to an input ROOT file
        isROOT (bool): True if input file is ROOT file, False if not (i.e. GDML)
    """
    foundTGeoManager, global_origins = ROOT.get_volume_position(input_file, isROOT, volName)
    if not foundTGeoManager:
        print(f'No TGeoManager found in {input_file}, cannot get module offsets.')
        return None
    elif foundTGeoManager:
        global_origins_list = [[global_origins.at(i).at(j) for j in range(global_origins.at(i).size())] for i in range(global_origins.size())]
        return global_origins_list

def main(input_file):
    if input_file.split('.')[-1] == 'root':
        isROOT = True
    else:
        isROOT = False
    ROOT.gROOT.ProcessLine('.L get_volume_position.cpp')
    # get positions of the all arapucas and place in csv file
    Range = [9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
    volumes = [f'volArapuca_0-{i}' for i in Range] + [f'volArapuca_1-{i}' for i in Range] \
		+ [f'volArapuca_2-{i}' for i in Range] + [f'volArapuca_3-{i}' for i in Range]
    with open('PDHD_PDS_Positions.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        for volName in volumes:
                module_offsets_GDML = get_volume_position(input_file, isROOT, volName)
                
                if module_offsets_GDML is not None and len(module_offsets_GDML) == 0:
                     print(f'Volume not found in TGeoManager of input file, check volume name in cpp file.')
                else:
                     print(f'volume position = {module_offsets_GDML}')
                     writer.writerow(module_offsets_GDML[0])
if __name__ == "__main__":
    fire.Fire(main)
