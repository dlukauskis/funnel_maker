import MDAnalysis as mda
import os
import numpy as np
import shutil
import subprocess as sp
import sys
import datetime

def name_file_sequentially(filename):
    """
    Checks if there is a file already with
    the same name. If yes, then gives the file
    a number in sequence. 
    
    Prevents file overwritting.
    
    Returns filename + sequential int (if needed) + file extension
    
    Eg. 
    
    1. os.listdir()
       ['file.txt']
    2. fname = name_file_sequentially('file.txt')
    3. print(fname)
       'file.0.txt'
    4. open(fname,'w').close
    5. os.listdir()
       ['file.txt','file.0.txt']
    
    Repeat 2-5
    ['file.txt','file.0.txt','file.1.txt']
    """
    
    if '/' not in filename:
        lst = [ f for f in os.listdir('.') if filename.split('.')[0] in f]
    else:
        path = os.path.abspath(filename)
        lst = [ f for f in os.listdir(path) if filename.split('.')[0] in f]
    
    if len(lst) == 0:
        
        return filename
    
    else:
        all_indices = [ int(f.split('.')[-2]) for f in lst if len(f.split('.')) > 2]
        if len(all_indices) == 0:
            new_last_index = 0
        else:
            last_index = max(all_indices)
            new_last_index = last_index + 1
        
    return filename.split('.')[0]+'.%i.'% (new_last_index)+filename.split('.')[-1]

def write_funnel_pymol_session(structure_file, p1, p2, topology_file='', 
                               extent = 0.85, extent_buffer = 0.15, 
                               l_proj = 0.5, u_proj = 3.5, beta_cent = 1.5, 
                               s_cent = 2.0, display_grid = False, center_grid_point=[], 
                               grid_display_lines=[], special_point=[]):
    
    """
    Writes a plumed script that loads the structure (and topology) file
    and draws a funnel based on your inputs.
    
    If provided with the data, can also display the search grid 
    used to define the funnel.
    """
    
    session_file_name = name_file_sequentially('visualise_funnel.pml')
    
    working_dir = os.path.dirname(structure_file)
    # some exceptions for amber files
    if 'rst' in structure_file.split('/')[-1]:
        top_prefix = topology_file.split('/')[-1].split('.')[0]
        new_topology_file = os.path.join(working_dir, top_prefix + '.top')
        
        strc_prefix = structure_file.split('/')[-1].split('.')[0]
        new_structure_file = os.path.join(working_dir, top_prefix + '.rst')
        
        shutil.copy(structure_file, new_structure_file)
        shutil.copy(topology_file, new_topology_file)
        
    p1_string = ''
    p2_string = ''
    
    for i in p1:
        p1_string += 'ID %i '%i
    
    for i in p2:
        p2_string += 'ID %i '%i
        
    with open(session_file_name, 'w') as PML_SCRIPT:
        if 'rst7' in structure_file.split('/')[-1]:
            PML_SCRIPT.write('load %s,mol\n'% new_topology_file)
            PML_SCRIPT.write('load %s,mol\n'% new_structure_file)
        else:
            PML_SCRIPT.write('load %s\n'% structure_file)
        PML_SCRIPT.write('\n')
        PML_SCRIPT.write('hide everythin, solvent\n')
        PML_SCRIPT.write('show cartoon, all\n')
        PML_SCRIPT.write('show sticks, organic\n')
        PML_SCRIPT.write('\n')
        n = 0
        if display_grid is True:    
            PML_SCRIPT.write('pseudoatom center_grid, pos=[%.3f, %.3f, %.3f]\n'%(center_grid_point[0],\
                                                                                 center_grid_point[1],\
                                                                                 center_grid_point[2]))
            PML_SCRIPT.write('show spheres, center_grid\n')
            PML_SCRIPT.write('pseudoatom special, pos=[%.3f, %.3f, %.3f]\n'%(special_point[0],\
                                                                             special_point[1],\
                                                                             special_point[2]))
            PML_SCRIPT.write('show spheres, special\n')
            PML_SCRIPT.write('color cyan, special\n')    
            for line in grid_display_lines:
                PML_SCRIPT.write(line)
        PML_SCRIPT.write('# take a point 10A deeper into the protein core\n')        
        PML_SCRIPT.write('select p1, %s\n'%p1_string)
        PML_SCRIPT.write('com p1\n')
        PML_SCRIPT.write('show spheres, p1_COM\n')
        PML_SCRIPT.write('\n')
        PML_SCRIPT.write('# select all CA within 10A of the grid center, find its COM\n')
        PML_SCRIPT.write('select p2, %s\n'%p2_string)
        PML_SCRIPT.write('com p2\n')
        PML_SCRIPT.write('show spheres, p2_COM\n')
        PML_SCRIPT.write('\n')
        PML_SCRIPT.write('\n')
        PML_SCRIPT.write('# draw the funnel\n')
        PML_SCRIPT.write('draw_funnel p1_COM, p2_COM, upper_wall=%f,\
                                                      lower_wall=%f,\
                                                      s_cent=%f,\
                                                      beta_cent=%f,\
                                                      wall_width=%f,\
                                                      wall_buffer=%f\n'% (u_proj*10,l_proj*10,s_cent*10, 
                                                                           beta_cent,extent*10,extent_buffer*10)) 
                                                     # x10 because its in A, not nm

def make_funnel(structure_file, topology_file='',
                grid_coords = [],
                ligand_name = 'MOL',
                output_pymol_session = False,
                display_grid = False):
    """
    Defines the funnel based on search criteria, essentially pointing away 
    from the protein.
    
    Inputs:
    
    structure_file - any MDAnalysis readable structure file
    topology_file - if using Amber structure files, you must provide a topology file as well
    grid_coords - if choosing to use xyz to define the centre of the binding site
                  (ala docking), provide the coords in a list
    ligand_name - if not providing the centre of the search grid in xyz, 
                    specify the name of the ligand. Default - MOL
    write_pymol_session - writes a pymol session, showing how the funnel was drawn
    display_grid - if outputing a pymol session, can project the search grid too. 
                   Basically for debugging.
                   
    Returns:
    
    p1 and p2 - numpy arrays that contain atom IDs as found in the structure file
                   
    """
    if '.rst7' in structure_file.split('/')[-1]:
        u = mda.Universe(topology_file,structure_file,format='RESTRT')
    elif '.inpcrd' in structure_file.split('/')[-1]:
        u = mda.Universe(topology_file,structure_file,format='INPCRD')
    else:
        u = mda.Universe(structure_file)
    
    if grid_coords:
        center_grid_point = grid_coords
    else:
        center_grid_point = u.select_atoms('resname %s'% ligand_name).center_of_mass()
        
        if len(center_grid_point) == 0:
            center_grid_point = u.select_atoms('resname LIG').center_of_mass()
            if len(center_grid_point) == 0:
                sys.exit("""Neither rename 'MOL' nor 'LIG' were found in the structure file.\n
                         Please specify either coordinates for the grid or an appropriate residue\n
                         name for the ligand.""")

    # lets project a grid on the binding site
    
    # 20 angstrom grid edge
    grid_length = 20
    # how many grid points along each edge
    n_edge_points = 5 
    # distance from each grid point to look for protein
    search_radius = (grid_length / n_edge_points) / 2 

    n = 0

    grid_display_lines = []
    solvent_com = []

    x_min, x_max = center_grid_point[0] - grid_length/2, center_grid_point[0] + grid_length/2
    y_min, y_max = center_grid_point[1] - grid_length/2, center_grid_point[1] + grid_length/2
    z_min, z_max = center_grid_point[2] - grid_length/2, center_grid_point[2] + grid_length/2 

    for x in np.linspace(x_min, x_max, n_edge_points):
        for y in np.linspace(y_min, y_max, n_edge_points):
            for z in np.linspace(z_min, z_max, n_edge_points):
                if display_grid is True:
                    grid_display_lines.append('pseudoatom g_pnt%i, pos=[%.3f, %.3f, %.3f]\n'%(n, x, y, z))
                    grid_display_lines.append('show spheres, g_pnt%i\n'%(n))

                protein_atoms = u.select_atoms('protein and point %.3f %.3f %.3f %.3f'% (x,y,z,\
                                                                                        search_radius),\
                                               periodic = False)
                if len(protein_atoms) > 0:
                    if display_grid is True:
                        grid_display_lines.append('color red, g_pnt%i\n'%(n))
                else:
                    solvent_com.append([x,y,z])
                n += 1            
    
    # if no protein detected at multiple grid points, find where the average coord is
    special_point = np.array([np.mean(np.array(solvent_com)[:,0]),\
                              np.mean(np.array(solvent_com)[:,1]),\
                              np.mean(np.array(solvent_com)[:,2])])

    # select all C alpha atoms within 10 A of the ligand/center of the grid
    near_ca = u.select_atoms('name CA and point %.3f %.3f %.3f 10'% (center_grid_point[0],\
                                                                    center_grid_point[1],\
                                                                    center_grid_point[2]),
                             periodic = False)

    initial_funnel_vector = special_point - near_ca.center_of_mass()
    initial_normed_funnel_vector = initial_funnel_vector / np.sqrt(np.sum(np.square(initial_funnel_vector)))

    into_the_protein = near_ca.center_of_mass() - 10 * initial_normed_funnel_vector

    p1_group = u.select_atoms('name CA and point %.3f %.3f %.3f 7'%(into_the_protein[0],\
                                                                    into_the_protein[1],\
                                                                    into_the_protein[2]),\
                             periodic = False)
    p1 = []
    p2 = []

    for i in p1_group.ids:
        p1.append(i)

    p1 = np.array(p1)

    for i in near_ca.ids:
        p2.append(i)

    p2 = np.array(p2)
    
    # build the 'real' funnel vector

    funnel_vector = near_ca.center_of_mass() - p1_group.center_of_mass()
    normed_funnel_vector = funnel_vector / np.sqrt(np.sum(np.square(funnel_vector)))

    funnel_vector_pnts = []

    if output_pymol_session is True:
        # this bit is for display purposes only...
        write_funnel_pymol_session(structure_file=structure_file, p1=p1, p2=p2, 
                                   topology_file=topology_file, display_grid=display_grid,
                                   center_grid_point=center_grid_point, 
                                   grid_display_lines=grid_display_lines, 
                                   special_point=special_point)
        
    return p1, p2
    
def get_protein_ligand_ids(structure_file, topology_file='', ligand_name = 'MOL'):
    """
    Returns a np array of atom IDs for 
    
    the protein and the ligand.
    
    Used for writing plumed.dat files
    """
    
    if '.rst7' in structure_file.split('/')[-1]:
        u = mda.Universe(topology_file,structure_file,format='RESTRT')
    elif '.inpcrd' in structure_file.split('/')[-1]:
        u = mda.Universe(topology_file,structure_file,format='INPCRD')
    else:
        u = mda.Universe(structure_file)
        
    ligand_group = u.select_atoms('resname %s'% ligand_name)
        
    if len(ligand_group) == 0:
        ligand_group = u.select_atoms('resname LIG')
        if len(ligand_group) == 0:
            sys.exit("""Neither rename 'MOL' nor 'LIG' were found in the structure file.\n
                      Please specify an appropriate residue name for the ligand.""")
    
    protein_group = u.select_atoms('protein')
    
    protein_IDs = []
    for i in protein_group.ids:
        protein_IDs.append(i)
    protein_IDs = np.array(protein_IDs)
    
    ligand_IDs = []
    for i in ligand_group.ids:
        ligand_IDs.append(i)
    ligand_IDs = np.array(ligand_IDs)
    
    return protein_IDs, ligand_IDs
    
def write_plumed_file(p1, p2, protein_IDs, lig_IDs, extent = 0.60, extent_buffer = 0.15, 
                      l_proj = 0.5, u_proj = 4.0, beta_cent = 1.5, 
                      s_cent = 2, deposition_pace = 1000,
                      print_pace = 1000, write_ProjectionOnAxis = False):
    
    """
    Writes a standard wt fun-metaD plumed.dat file in the current working directory.
    
    p1, p2 - numpy array, atom IDs that will act as anchor points for the funnel
    
    protein_IDs - numpy array, atom IDs (inclusive) belonging to the protein / host molecule
    lig_IDs - numpy array, atom IDs (inclusive) belonging to the ligand / guest molecule
    
    Length units are in nm.
    """
    version = 1.0
    
    p1_str = ''
    for i in p1:
        p1_str += str(i) + ','
    
    p1_str = p1_str[:-1]

    p2_str = ''
    for i in p2:
        p2_str += str(i) + ','
    
    p2_str = p2_str[:-1]
    
    protein_str = '%i-%i'% (protein_IDs[0], protein_IDs[-1])
    lig_str = '%i-%i'% (lig_IDs[0], lig_IDs[-1])

    with open('plumed.dat', 'w') as FILE:
        FILE.write('####################################\n')
        FILE.write('#plumed.dat for Funnel Metadynamics#\n')
        FILE.write('# Written on %s\n'% datetime.datetime.now())
        FILE.write('# By funnel_maker %s\n'% str(version))
        FILE.write('####################################\n')
        FILE.write('RESTART\n')
        FILE.write('\n')
        FILE.write('###############################################\n')
        FILE.write('###DEFINE RADIUS + CALC PROT-LIG VECTOR COMP###\n')
        FILE.write('###############################################\n')
        if write_ProjectionOnAxis is True:
            FILE.write('LOAD FILE=ProjectionOnAxis.cpp\n')
        FILE.write('\n')
        FILE.write('WHOLEMOLECULES STRIDE=1 ENTITY0=%s ENTITY1=%s\n'% (protein_str, lig_str))
        FILE.write('\n')
        FILE.write('########################\n')
        FILE.write('###DEFINITION_OF_COMs###\n')
        FILE.write('########################\n')
        FILE.write('lig: COM ATOMS=%s\n'% lig_str)
        FILE.write('p1: COM ATOMS=%s\n'% p1_str)
        FILE.write('p2: COM ATOMS=%s\n'% p2_str)
        FILE.write('\n')
        FILE.write('\n')
        FILE.write('########################\n')
        FILE.write('###DEFINITION_OF_ARGs###\n')
        FILE.write('########################\n')
        FILE.write('# CV1: pp.proj = projection on the axis. The distance from the axis to the origin (along the axis)\n')
        FILE.write('# CV2: pp.ext = orthogonal distance between the ATOM(=lig) to the axis.\n')
        FILE.write('\n')
        FILE.write('############\n')
        FILE.write('###PoA_CV ##\n')
        FILE.write('############\n')
        FILE.write('pp: PROJECTION_ON_AXIS AXIS_ATOMS=p1,p2 ATOM=lig\n')
        FILE.write('\n')
        FILE.write('#######################\n')
        FILE.write('###FUNNEL_PARAMETERS###\n')
        FILE.write('#######################\n')
        FILE.write('s_cent: CONSTANT VALUES=%.1f                                    # INFLEXION\n'% s_cent)
        FILE.write('beta_cent: CONSTANT VALUES=%.1f                                 # STEEPNESS\n'% beta_cent)
        FILE.write('wall_width: CONSTANT VALUES=%.2f                                # WIDTH (h)\n'% extent)
        FILE.write('wall_buffer: CONSTANT VALUES=%.2f                               # BUFFER (f, total width = WIDTH + BUFFER)\n'% extent_buffer)
        FILE.write('lwall: LOWER_WALLS ARG=pp.proj AT=%.1f KAPPA=20000.0 EXP=2 EPS=1 # Lower Wall (the starting point of the funnel)\n'% l_proj)
        FILE.write('uwall: UPPER_WALLS ARG=pp.proj AT=%.1f KAPPA=20000.0 EXP=2 EPS=1 # Upper Wall (the ending point of the funnel)\n'% u_proj)
        FILE.write('\n')
        FILE.write('##################################\n')
        FILE.write('###########CALCULATE FUNNEL#######\n')
        FILE.write('# Returns the radius of the funnel\n')
        FILE.write('# at the current value of the cv\n')
        FILE.write('##################################\n')
        FILE.write('MATHEVAL ...\n')
        FILE.write('        LABEL=wall_center\n')
        FILE.write('        ARG=pp.proj,s_cent,beta_cent,wall_width,wall_buffer\n')
        FILE.write('        VAR=s,sc,b,h,f\n')
        FILE.write('        FUNC=h*(1./(1.+exp(b*(s-sc))))+f\n')
        FILE.write('        PERIODIC=NO\n')
        FILE.write('... MATHEVAL\n')
        FILE.write('\n')
        FILE.write('\n')
        FILE.write('##############################\n')
        FILE.write('#####POTENTIAL_PARAMETERS#####\n')
        FILE.write('##############################\n')
        FILE.write('scaling: CONSTANT VALUES=1.0\n')
        FILE.write('spring: CONSTANT VALUES=1000.0\n')
        FILE.write('\n')
        FILE.write('##############################\n')
        FILE.write('#######DEFINE_POTENTIAL#######\n')
        FILE.write('##############################\n')
        FILE.write('MATHEVAL ...\n')
        FILE.write('        LABEL=wall_bias\n')
        FILE.write('        ARG=pp.ext,spring,wall_center,scaling\n')
        FILE.write('        VAR=z,k,zc,sf\n')
        FILE.write('        FUNC=step(z-zc)*k*(z-zc)*(z-zc)/(sf*sf)\n')
        FILE.write('        PERIODIC=NO\n')
        FILE.write('... MATHEVAL\n')
        FILE.write('\n')
        FILE.write('finalbias: BIASVALUE ARG=wall_bias\n')
        FILE.write('\n')
        FILE.write('\n')
        FILE.write('###############################\n')
        FILE.write('########DEFINE_METAD###########\n')
        FILE.write('###############################\n')
        FILE.write('METAD ...\n')
        FILE.write('        LABEL=meta ARG=pp.proj,pp.ext \n')
        FILE.write('        SIGMA=0.025,0.03 HEIGHT=1.5 \n')
        FILE.write('        PACE=%i FILE=HILLS \n'% deposition_pace)
        FILE.write('        GRID_MIN=%.1f,%.1f GRID_MAX=%.1f,%.1f GRID_SPACING=0.005,0.06\n'% \
                                       ((l_proj - 0.5),0.0, # proj min, extent min
                                        (u_proj + 0.5),(extent+extent_buffer+0.2))) # proj max, extent max
        FILE.write('        BIASFACTOR=10.0 TEMP=298\n')
        FILE.write('... METAD\n')
        FILE.write('\n')
        FILE.write('PRINT ARG=* STRIDE=%i FILE=COLVAR FMT=%%8.4f\n'% print_pace)
        
def write_opes_plumed_file(p1, p2, protein_IDs, lig_IDs, extent = 0.60, extent_buffer = 0.15, 
                      l_proj = 0.5, u_proj = 4.0, beta_cent = 1.5, 
                      s_cent = 2, deposition_pace = 1000, barrier = 60,
                      print_pace = 1000, write_ProjectionOnAxis = False):
    
    """
    Writes a standard wt fun-OPES plumed.dat file
    
    p1, p2 - numpy array, atom IDs that will act as anchor points for the funnel
    
    protein_IDs - numpy array, atom IDs (inclusive) belonging to the protein / host molecule
    lig_IDs - numpy array, atom IDs (inclusive) belonging to the ligand / guest molecule
    
    Length units are in nm.
    """
    version = 1.0
    
    p1_str = ''
    for i in p1:
        p1_str += str(i) + ','
    
    p1_str = p1_str[:-1]

    p2_str = ''
    for i in p2:
        p2_str += str(i) + ','
    
    p2_str = p2_str[:-1]
    
    protein_str = '%i-%i'% (protein_IDs[0], protein_IDs[-1])
    lig_str = '%i-%i'% (lig_IDs[0], lig_IDs[-1])
    
    with open('plumed.dat', 'w') as FILE:
        FILE.write('####################################\n')
        FILE.write('#plumed.dat for Funnel OPES#\n')
        FILE.write('# Written on %s\n'% datetime.datetime.now())
        FILE.write('# By funnel_maker %s\n'% str(version))
        FILE.write('####################################\n')
        FILE.write('RESTART\n')
        FILE.write('\n')
        FILE.write('###############################################\n')
        FILE.write('###DEFINE RADIUS + CALC PROT-LIG VECTOR COMP###\n')
        FILE.write('###############################################\n')
        if write_ProjectionOnAxis is True:
            FILE.write('LOAD FILE=ProjectionOnAxis.cpp\n')
        FILE.write('LOAD FILE=OPESwt.cpp\n')
        FILE.write('\n')
        FILE.write('WHOLEMOLECULES STRIDE=1 ENTITY0=%s ENTITY1=%s\n'% (protein_str, lig_str))
        FILE.write('\n')
        FILE.write('########################\n')
        FILE.write('###DEFINITION_OF_COMs###\n')
        FILE.write('########################\n')
        FILE.write('lig: COM ATOMS=%s\n'% lig_str)
        FILE.write('p1: COM ATOMS=%s\n'% p1_str)
        FILE.write('p2: COM ATOMS=%s\n'% p2_str)
        FILE.write('\n')
        FILE.write('\n')
        FILE.write('########################\n')
        FILE.write('###DEFINITION_OF_ARGs###\n')
        FILE.write('########################\n')
        FILE.write('# CV1: pp.proj = projection on the axis. The distance from the axis to the origin (along the axis)\n')
        FILE.write('# CV2: pp.ext = orthogonal distance between the ATOM(=lig) to the axis.\n')
        FILE.write('\n')
        FILE.write('############\n')
        FILE.write('###PoA_CV ##\n')
        FILE.write('############\n')
        FILE.write('pp: PROJECTION_ON_AXIS AXIS_ATOMS=p1,p2 ATOM=lig\n')
        FILE.write('\n')
        FILE.write('#######################\n')
        FILE.write('###FUNNEL_PARAMETERS###\n')
        FILE.write('#######################\n')
        FILE.write('s_cent: CONSTANT VALUES=%.1f                                    # INFLEXION\n'% s_cent)
        FILE.write('beta_cent: CONSTANT VALUES=%.1f                                 # STEEPNESS\n'% beta_cent)
        FILE.write('wall_width: CONSTANT VALUES=%.2f                                # WIDTH (h)\n'% extent)
        FILE.write('wall_buffer: CONSTANT VALUES=%.2f                               # BUFFER (f, total width = WIDTH + BUFFER)\n'% extent_buffer)
        FILE.write('lwall: LOWER_WALLS ARG=pp.proj AT=%.1f KAPPA=20000.0 EXP=2 EPS=1 # Lower Wall (the starting point of the funnel)\n'% l_proj)
        FILE.write('uwall: UPPER_WALLS ARG=pp.proj AT=%.1f KAPPA=20000.0 EXP=2 EPS=1 # Upper Wall (the ending point of the funnel)\n'% u_proj)
        FILE.write('\n')
        FILE.write('##################################\n')
        FILE.write('###########CALCULATE FUNNEL#######\n')
        FILE.write('# Returns the radius of the funnel\n')
        FILE.write('# at the current value of the cv\n')
        FILE.write('##################################\n')
        FILE.write('MATHEVAL ...\n')
        FILE.write('        LABEL=wall_center\n')
        FILE.write('        ARG=pp.proj,s_cent,beta_cent,wall_width,wall_buffer\n')
        FILE.write('        VAR=s,sc,b,h,f\n')
        FILE.write('        FUNC=h*(1./(1.+exp(b*(s-sc))))+f\n')
        FILE.write('        PERIODIC=NO\n')
        FILE.write('... MATHEVAL\n')
        FILE.write('\n')
        FILE.write('\n')
        FILE.write('##############################\n')
        FILE.write('#####POTENTIAL_PARAMETERS#####\n')
        FILE.write('##############################\n')
        FILE.write('scaling: CONSTANT VALUES=1.0\n')
        FILE.write('spring: CONSTANT VALUES=1000.0\n')
        FILE.write('\n')
        FILE.write('##############################\n')
        FILE.write('#######DEFINE_POTENTIAL#######\n')
        FILE.write('##############################\n')
        FILE.write('MATHEVAL ...\n')
        FILE.write('        LABEL=wall_bias\n')
        FILE.write('        ARG=pp.ext,spring,wall_center,scaling\n')
        FILE.write('        VAR=z,k,zc,sf\n')
        FILE.write('        FUNC=step(z-zc)*k*(z-zc)*(z-zc)/(sf*sf)\n')
        FILE.write('        PERIODIC=NO\n')
        FILE.write('... MATHEVAL\n')
        FILE.write('\n')
        FILE.write('finalbias: BIASVALUE ARG=wall_bias\n')
        FILE.write('\n')
        FILE.write('\n')
        FILE.write('###############################\n')
        FILE.write('########DEFINE_OPES_WT#########\n')
        FILE.write('###############################\n')
        FILE.write('OPES_WT ...\n')
        FILE.write('  LABEL=opes\n')
        FILE.write('  FILE=Kernels.data\n')
        FILE.write('  TEMP=298\n')
        FILE.write('  ARG=pp.proj,pp.ext\n')
        FILE.write('  PACE=%s'% deposition_pace)
        FILE.write('  BARRIER=%i\n'% barrier)
        FILE.write('  PROB_WFILE=Prob.data\n')
        FILE.write('  PROB_WSTRIDE=50000\n')
        FILE.write('... OPES_WT\n')
        FILE.write('\n')
        FILE.write('PRINT ARG=* STRIDE=%i FILE=COLVAR FMT=%%8.4f\n'% print_pace)
        
def get_funnel_definitions_from_plumed(plumed_file):
    
    """
    Reads a 'plumed.dat' file and extracts the funnel definitions
    
    Returns:
    
    p1 - numpy array, atom IDs
    p2 - numpy array, atoms IDs
    s_cent - float
    beta_cent - float
    wall_width - float
    wall_buffer - float
    lwall - float
    uwall - float
    """
    
    plumed_file_lst = [line[:-1] for line in open(plumed_file,'r').readlines()]

    # s_cent: CONSTANT VALUES=2.1                                           # INFLEXION
    # beta_cent: CONSTANT VALUES=1.5    	                              # STEEPNESS
    # wall_width: CONSTANT VALUES=1.55                                      # WIDTH (h)
    # wall_buffer: CONSTANT VALUES=0.15  	                              # BUFFER (f, total width = WIDTH + BUFFER)
    # lwall: LOWER_WALLS ARG=pp.proj AT=0.5 KAPPA=20000.0 EXP=2 EPS=1        # Lower Wall (the starting point of the funnel)
    # uwall: UPPER_WALLS ARG=pp.proj AT=3.0 KAPPA=20000.0 EXP=2 EPS=1        # Upper Wall (the ending point of the funnel)

    
    p1 = []
    p2 = []
    for line in plumed_file_lst:
        if 'p1:' in line:
            p1_str = line.split('=')[1]

            p1_list = p1_str.split(',')
            for i in p1_list:
                p1.append(int(i))
                
            p1 = np.array(p1)
        elif 'p2:' in line:
            p2_str = line.split('=')[1]
            p2_list = p2_str.split(',')
            for i in p2_list:
                p2.append(int(i))
                
            p2 = np.array(p2)
        elif 's_cent:' in line:
            s_cent = float(line.split('=')[1][:4])
        elif 'beta_cent:' in line:
            beta_cent = float(line.split('=')[1][:4])
        elif 'wall_width:' in line:
            wall_width = float(line.split('=')[1][:4])
        elif 'wall_buffer:' in line:
            wall_buffer = float(line.split('=')[1][:4])
        elif 'lwall:' in line:
            lwall = float(line.split('=')[2][:3])
        elif 'uwall:' in line:
            uwall = float(line.split('=')[2][:3])

    return p1, p2, s_cent, beta_cent, wall_width, wall_buffer, lwall, uwall

def write_ProjectionOnAxis_script():

    with open('ProjectionOnAxis.cpp','w') as SCRIPT:
        SCRIPT.write("/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
        SCRIPT.write("    Copyright (c) 2011-2018 The plumed team\n")
        SCRIPT.write("    (see the PEOPLE file at the root of the distribution for a list of names)\n")
        SCRIPT.write("\n")
        SCRIPT.write("    See http://www.plumed.org for more information.\n")
        SCRIPT.write("\n")
        SCRIPT.write("    This file is part of plumed, version 2.\n")
        SCRIPT.write("\n")
        SCRIPT.write("    plumed is free software: you can redistribute it and/or modify\n")
        SCRIPT.write("    it under the terms of the GNU Lesser General Public License as published by\n")
        SCRIPT.write("    the Free Software Foundation, either version 3 of the License, or\n")
        SCRIPT.write("    (at your option) any later version.\n")
        SCRIPT.write("\n")
        SCRIPT.write("    plumed is distributed in the hope that it will be useful,\n")
        SCRIPT.write("    but WITHOUT ANY WARRANTY; without even the implied warranty of\n")
        SCRIPT.write("    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n")
        SCRIPT.write("    GNU Lesser General Public License for more details.\n")
        SCRIPT.write("\n")
        SCRIPT.write("    You should have received a copy of the GNU Lesser General Public License\n")
        SCRIPT.write("    along with plumed.  If not, see <http://www.gnu.org/licenses/>.\n")
        SCRIPT.write("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */\n")
        SCRIPT.write('#include "colvar/Colvar.h"\n')
        SCRIPT.write('#include "core/ActionRegister.h"\n')
        SCRIPT.write('#include "tools/Angle.h"\n')
        SCRIPT.write("\n")
        SCRIPT.write("#include <string>\n")
        SCRIPT.write("#include <cmath>\n")
        SCRIPT.write("\n")
        SCRIPT.write("using namespace std;\n")
        SCRIPT.write("\n")
        SCRIPT.write("namespace PLMD {\n")
        SCRIPT.write("namespace colvar {\n")
        SCRIPT.write("\n")
        SCRIPT.write("//+PLUMEDOC COLVAR PROJECTION_ON_AXIS\n")
        SCRIPT.write("/*\n")
        SCRIPT.write("Calculate the projection on an axis.\n")
        SCRIPT.write("\n")
        SCRIPT.write("\par Examples\n")
        SCRIPT.write("\n")
        SCRIPT.write("*/\n")
        SCRIPT.write("//+ENDPLUMEDOC\n")
        SCRIPT.write("\n")
        SCRIPT.write("class ProjectionOnAxis : public Colvar {\n")
        SCRIPT.write("  bool pbc;\n")
        SCRIPT.write("\n")
        SCRIPT.write("public:\n")
        SCRIPT.write("  explicit ProjectionOnAxis(const ActionOptions&);\n")
        SCRIPT.write("// active methods:\n")
        SCRIPT.write("  virtual void calculate();\n")
        SCRIPT.write("  static void registerKeywords( Keywords& keys );\n")
        SCRIPT.write("};\n")
        SCRIPT.write("\n")
        SCRIPT.write('PLUMED_REGISTER_ACTION(ProjectionOnAxis,"PROJECTION_ON_AXIS")\n')
        SCRIPT.write("\n")
        SCRIPT.write("void ProjectionOnAxis::registerKeywords( Keywords& keys ) {\n")
        SCRIPT.write("  Colvar::registerKeywords(keys);\n")
        SCRIPT.write('  keys.add("atoms","AXIS_ATOMS","the atoms that define the direction of the axis of interest");\n')
        SCRIPT.write('  keys.add("atoms","ATOM","the atom whose position we want to project on the axis of interest");\n')
        SCRIPT.write("}\n")
        SCRIPT.write("\n")
        SCRIPT.write("ProjectionOnAxis::ProjectionOnAxis(const ActionOptions&ao):\n")
        SCRIPT.write("  PLUMED_COLVAR_INIT(ao),\n")
        SCRIPT.write("  pbc(true)\n")
        SCRIPT.write("{\n")
        SCRIPT.write("  vector<AtomNumber> axis_atoms;\n")
        SCRIPT.write('  parseAtomList("AXIS_ATOMS",axis_atoms);\n')
        SCRIPT.write('  if( axis_atoms.size()!=2 ) error("should only be two atoms specified to AXIS_ATOMS keyword");\n')
        SCRIPT.write("  vector<AtomNumber> atom; \n")
        SCRIPT.write('  parseAtomList("ATOM",atom);\n')
        SCRIPT.write('  if( atom.size()!=1 ) error("should only be one atom specified to ATOM keyword");\n')
        SCRIPT.write('  log.printf("  calculating projection of vector connecting atom %d and atom %d on vector connecting atom %d and atom %d \\n",\n')
        SCRIPT.write("                axis_atoms[0].serial(), atom[0].serial(), axis_atoms[0].serial(), axis_atoms[1].serial() ); \n")
        SCRIPT.write("  bool nopbc=!pbc;\n")
        SCRIPT.write('  parseFlag("NOPBC",nopbc);\n')
        SCRIPT.write("  pbc=!nopbc;\n")
        SCRIPT.write("\n")
        SCRIPT.write('  if(pbc) log.printf("  using periodic boundary conditions\\n");\n')
        SCRIPT.write('  else    log.printf("  not using periodic boundary conditions\\n");\n')
        SCRIPT.write("\n")
        SCRIPT.write("  // Add values to store data\n")
        SCRIPT.write('  addComponentWithDerivatives("proj"); componentIsNotPeriodic("proj"); \n')
        SCRIPT.write('  addComponentWithDerivatives("ext"); componentIsNotPeriodic("ext");\n')
        SCRIPT.write("  // Get all the atom positions \n")
        SCRIPT.write("  axis_atoms.push_back( atom[0] ); \n")
        SCRIPT.write("  requestAtoms(axis_atoms);\n")
        SCRIPT.write("  checkRead();\n")
        SCRIPT.write("}\n")
        SCRIPT.write("// calculator\n")
        SCRIPT.write("void ProjectionOnAxis::calculate() {\n")
        SCRIPT.write("\n")
        SCRIPT.write("  Vector rik, rjk;\n")
        SCRIPT.write("  if( pbc ) {\n")
        SCRIPT.write("      rik = pbcDistance( getPosition(2), getPosition(0) );\n")
        SCRIPT.write("      rjk = pbcDistance( getPosition(2), getPosition(1) );\n")
        SCRIPT.write("  } else {\n")
        SCRIPT.write("      rik = delta( getPosition(2), getPosition(0) );\n")
        SCRIPT.write("      rjk = delta( getPosition(2), getPosition(1) );\n")
        SCRIPT.write("  }\n")
        SCRIPT.write("  Vector rij = delta( rik, rjk ); double dij = rij.modulo(); \n")
        SCRIPT.write("  Vector nij = (1.0/dij)*rij; Tensor dij_a1;\n")
        SCRIPT.write("  // Derivative of director connecting atom1 - atom2 wrt the position of atom 1\n")
        SCRIPT.write("  dij_a1(0,0) = ( -(nij[1]*nij[1]+nij[2]*nij[2])/dij );   // dx/dx\n")
        SCRIPT.write("  dij_a1(0,1) = (  nij[0]*nij[1]/dij );                   // dx/dy\n")
        SCRIPT.write("  dij_a1(0,2) = (  nij[0]*nij[2]/dij );                   // dx/dz\n")
        SCRIPT.write("  dij_a1(1,0) = (  nij[1]*nij[0]/dij );                   // dy/dx\n")
        SCRIPT.write("  dij_a1(1,1) = ( -(nij[0]*nij[0]+nij[2]*nij[2])/dij );   // dy/dy\n")
        SCRIPT.write("  dij_a1(1,2) = (  nij[1]*nij[2]/dij );\n")
        SCRIPT.write("  dij_a1(2,0) = (  nij[2]*nij[0]/dij );\n")
        SCRIPT.write("  dij_a1(2,1) = (  nij[2]*nij[1]/dij );\n")
        SCRIPT.write("  dij_a1(2,2) = ( -(nij[1]*nij[1]+nij[0]*nij[0])/dij );\n")
        SCRIPT.write("\n")
        SCRIPT.write("  // Calculate dot product and derivatives\n")
        SCRIPT.write("  double d = dotProduct( -rik, nij );\n")
        SCRIPT.write("  Vector dd1 = matmul(-rik, dij_a1) - nij;\n")
        SCRIPT.write("  Vector dd2 = matmul(rik, dij_a1);\n")
        SCRIPT.write("  Vector dd3 = nij;\n")
        SCRIPT.write('  Value* pval=getPntrToComponent("proj"); pval->set( d );\n')
        SCRIPT.write("  setAtomsDerivatives( pval, 0, dd1 );\n")
        SCRIPT.write("  setAtomsDerivatives( pval, 1, dd2 );\n")
        SCRIPT.write("  setAtomsDerivatives( pval, 2, dd3 );\n")
        SCRIPT.write("  setBoxDerivatives( pval, -Tensor( rik, dd1 ) - Tensor( rjk, dd2 ) );\n")
        SCRIPT.write("  // Calculate derivatives of perpendicular distance from axis\n")
        SCRIPT.write("  double c = sqrt( rik.modulo2() - d*d ); double invc = (1.0/c);\n")
        SCRIPT.write("  // Calculate derivatives of the other thing\n")
        SCRIPT.write("  Vector der1 = invc*(rik - d*dd1);\n")
        SCRIPT.write("  Vector der2 = invc*(-d*dd2);\n")
        SCRIPT.write("  Vector der3 = invc*(-rik - d*dd3);\n")
        SCRIPT.write("\n")
        SCRIPT.write('  Value* cval=getPntrToComponent("ext"); cval->set( c ); \n')
        SCRIPT.write("  setAtomsDerivatives( cval, 0, der1 );\n")
        SCRIPT.write("  setAtomsDerivatives( cval, 1, der2 );\n")
        SCRIPT.write("  setAtomsDerivatives( cval, 2, der3 );\n")
        SCRIPT.write("  setBoxDerivatives( cval, -Tensor( rik, der1 ) - Tensor( rjk, der2 ) );\n")
        SCRIPT.write("}\n")
        SCRIPT.write("\n")
        SCRIPT.write("}\n")
        SCRIPT.write("}\n")

def write_run_py(structure_file,topology_file,run_time,lig_ids,p1_ids,p2_ids,lower_wall,upper_wall,wall_buffer,wall_width,s_cent,beta_cent):
    with open('run.py','w') as FILE:
        FILE.write('from simtk.openmm import *\n')
        FILE.write('from simtk.openmm.app import *\n')
        FILE.write('from simtk.unit import *\n')
        FILE.write('from simtk.openmm.app.metadynamics import *\n')
        FILE.write('from parmed import load_file\n')
        FILE.write('from sys import stdout\n')
        FILE.write('import os\n')
        FILE.write('import subprocess as sp\n')
        FILE.write('import shutil\n')
        FILE.write('\n')
        FILE.write('import numpy as np\n')
        FILE.write("import MDAnalysis as mda\n")
        FILE.write('\n')
        FILE.write('# if you want to assign which GPU to use,\n')
        FILE.write('# uncomment the line below\n')
        FILE.write("#os.environ['CUDA_VISIBLE_DEVICES'] = '0'\n")
        FILE.write('\n')
        FILE.write("\n")
        FILE.write('"""\n')
        FILE.write('OpenMM-native implementation of funnel-metadynamics\n')
        FILE.write('Steps:\n')
        FILE.write('    1. 10k enmin steps\n')
        FILE.write('    2. 5ns solute-restrained NPT equilibration, with MC bstat\n')
        FILE.write('    3. 5ns Calpha+ligand-restrained NVT equilibration\n')
        FILE.write(f'    4. {run_time}ns well-tempered funnel-metadynamics\n')
        FILE.write('\n')
        FILE.write('Outputs:\n')
        FILE.write('    1. COLVAR, logging CVs every 2 ps\n')
        FILE.write('    2. bias files, written every 1 ns\n')
        FILE.write('    3. trajectory, written every 100 ps\n')
        FILE.write('    4. checkpoint, every 100 ps\n')
        FILE.write('"""\n')
        FILE.write("\n")
        FILE.write("def get_CA_indices(coords, parm):\n")
        FILE.write("\n")
        FILE.write("    if coords.endswith('.gro'):\n")
        FILE.write("        u = mda.Universe(coords)\n")
        FILE.write("    elif coords.endswith('.rst') or coords.endswith('.rst7') or coords.endswith('.inpcrd'):\n")
        FILE.write("\n")
        FILE.write("        u = mda.Universe(parm, coords, format='INPCRD')\n")
        FILE.write("\n")
        FILE.write("    ca_atoms = u.select_atoms('name CA')\n")
        FILE.write("\n")
        FILE.write("    ca_list = [int(ca) for ca in ca_atoms.indices]\n")
        FILE.write("\n")
        FILE.write("    return ca_list\n")
        FILE.write("\n")
        FILE.write("def get_ligand_indices(coords, parm, lig_name = 'MOL', include_h = True):\n")
        FILE.write("\n")
        FILE.write("    if coords.endswith('.gro'):\n")
        FILE.write("        u = mda.Universe(coords)\n")
        FILE.write("    elif coords.endswith('.rst') or coords.endswith('.rst7') or coords.endswith('.inpcrd'):\n")
        FILE.write("\n")
        FILE.write("        u = mda.Universe(parm, coords, format='INPCRD')\n")
        FILE.write("\n")
        FILE.write("    if include_h is False:\n")
        FILE.write("        lig_atoms = u.select_atoms('resname %s and not (name h* or name H*)'% lig_name)\n")
        FILE.write("    else:\n")
        FILE.write("        lig_atoms = u.select_atoms('resname %s'% lig_name)\n")
        FILE.write("\n")
        FILE.write("    lig_list = [int(a) for a in lig_atoms.indices]\n")
        FILE.write("\n")
        FILE.write("    return lig_list\n")
        FILE.write("\n")
        FILE.write("def get_protein_indices(coords, parm, include_h = True):\n")
        FILE.write("\n")
        FILE.write("    if coords.endswith('.gro'):\n")
        FILE.write("        u = mda.Universe(coords)\n")
        FILE.write("    elif coords.endswith('.rst') or coords.endswith('.rst7') or coords.endswith('.inpcrd'):\n")
        FILE.write("\n")
        FILE.write("        u = mda.Universe(parm, coords, format='INPCRD')\n")
        FILE.write("\n")
        FILE.write("    if include_h is False:\n")
        FILE.write("        prot_atoms = u.select_atoms('protein and not name H*')\n")
        FILE.write("    else:\n")
        FILE.write("        prot_atoms = u.select_atoms('protein')\n")
        FILE.write("\n")
        FILE.write("    prot_list = [int(a) for a in prot_atoms.indices]\n")
        FILE.write("\n")
        FILE.write("    return prot_list\n")
        FILE.write("\n")
        FILE.write("def run_10k_enmin(params, input_positions):\n")
        FILE.write("\n")
        FILE.write("    # prepare system\n")
        FILE.write("    system = params.createSystem(nonbondedMethod=PME,\n")
        FILE.write("                               nonbondedCutoff=1.0*nanometers,\n")
        FILE.write("                               constraints=HBonds,\n")
        FILE.write("                               rigidWater=True,\n")
        FILE.write("                               ewaldErrorTolerance=0.0005)\n")
        FILE.write("\n")
        FILE.write("    integrator = LangevinIntegrator(300*kelvin, 1.0/picoseconds,\n")
        FILE.write("                                      2.0*femtoseconds)\n")
        FILE.write("    integrator.setConstraintTolerance(0.00001)\n")
        FILE.write("\n")
        FILE.write("    # prepare simulation\n")
        FILE.write("    platform = Platform.getPlatformByName('CUDA')\n")
        FILE.write("    properties = {'CudaPrecision':'mixed'}\n")
        FILE.write("    simulation = Simulation(params.topology, system,\n")
        FILE.write("                            integrator, platform, properties)\n")
        FILE.write("\n")
        FILE.write("    simulation.context.setPositions(input_positions)\n")
        FILE.write("\n")
        FILE.write("    ### 10k min, with 10kJ target ###\n")
        FILE.write("    simulation.minimizeEnergy(maxIterations=10000, tolerance=10*kilojoule/mole)\n")
        FILE.write("\n")
        FILE.write("    p = simulation.context.getState(getPositions=True).getPositions()\n")
        FILE.write("\n")
        FILE.write("    return p\n")
        FILE.write("\n")
        FILE.write("def run_5ns_NPT_restaint_equil(params, input_positions):\n")
        FILE.write("\n")
        FILE.write("    ### 5ns in NPT (MonteCarlo bstat) with 5kcal all heavy solute atom restraints ###\n")
        FILE.write("    # prepare system\n")
        FILE.write("    system = params.createSystem(nonbondedMethod=PME,\n")
        FILE.write("                               nonbondedCutoff=1.0*nanometers,\n")
        FILE.write("                               constraints=HBonds,\n")
        FILE.write("                               rigidWater=True,\n")
        FILE.write("                               ewaldErrorTolerance=0.0005)\n")
        FILE.write("\n")
        FILE.write("    restraint = HarmonicBondForce()\n")
        FILE.write("    restraint.setUsesPeriodicBoundaryConditions(True)\n")
        FILE.write("\n")
        FILE.write("    system.addForce(restraint)\n")
        FILE.write("    nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]\n")
        FILE.write("    atomsToRestrain = get_protein_indices(coords_file, parm_file, include_h = False) + get_ligand_indices(coords_file, parm_file, include_h = False)\n")
        FILE.write("\n")
        FILE.write("    dummyIndex = []\n")
        FILE.write("    positions = input_positions\n")
        FILE.write("    for i in atomsToRestrain:\n")
        FILE.write("        j = system.addParticle(0)\n")
        FILE.write("        nonbonded.addParticle(0, 1, 0)\n")
        FILE.write("        nonbonded.addException(i, j, 0, 1, 0)\n")
        FILE.write("        restraint.addBond(i, j, 0*nanometers, 5*kilocalories_per_mole/angstrom**2)\n")
        FILE.write("        dummyIndex.append(j)\n")
        FILE.write("        input_positions.append(positions[i])\n")
        FILE.write("    integrator = LangevinIntegrator(300*kelvin, 1.0/picoseconds,\n")
        FILE.write("                                      2.0*femtoseconds)\n")
        FILE.write("    system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))\n")
        FILE.write("    context = Context(system, integrator)\n")
        FILE.write("\n")
        FILE.write("    context.setPositions(input_positions)\n")
        FILE.write("\n")
        FILE.write("    run_time = 5.0 # ns\n")
        FILE.write("\n")
        FILE.write("    print('Initial energy:', context.getState(getEnergy=True).getPotentialEnergy())\n")
        FILE.write("    integrator.step(int(run_time * 500000))\n")
        FILE.write("    print('Final energy:', context.getState(getEnergy=True).getPotentialEnergy())\n")
        FILE.write("\n")
        FILE.write("    p = context.getState(getPositions=True).getPositions()[:dummyIndex[0]]\n")
        FILE.write("\n")
        FILE.write("    return p\n")
        FILE.write("\n")
        FILE.write("def run_5ns_NVT_restaint_equil(params, input_positions):\n")
        FILE.write("\n")
        FILE.write("    ### 5ns in NVT with 5kcal for Ca and ligand heavy atoms ###\n")
        FILE.write("\n")
        FILE.write("    # prepare system\n")
        FILE.write("    system = params.createSystem(nonbondedMethod=PME,\n")
        FILE.write("                                 nonbondedCutoff=1.0*nanometers,\n")
        FILE.write("                                 constraints=HBonds,\n")
        FILE.write("                                 rigidWater=True,\n")
        FILE.write("                                 ewaldErrorTolerance=0.0005)\n")
        FILE.write("\n")
        FILE.write("    restraint = HarmonicBondForce()\n")
        FILE.write("    restraint.setUsesPeriodicBoundaryConditions(True)\n")
        FILE.write("\n")
        FILE.write("    system.addForce(restraint)\n")
        FILE.write("    nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]\n")
        FILE.write("    atomsToRestrain = get_CA_indices(coords_file, parm_file) + get_ligand_indices(coords_file, parm_file, include_h = False)\n")
        FILE.write("\n")
        FILE.write("    dummyIndex = []\n")
        FILE.write("    positions = input_positions\n")
        FILE.write("    for i in atomsToRestrain:\n")
        FILE.write("        j = system.addParticle(0)\n")
        FILE.write("        nonbonded.addParticle(0, 1, 0)\n")
        FILE.write("        nonbonded.addException(i, j, 0, 1, 0)\n")
        FILE.write("        restraint.addBond(i, j, 0*nanometers, 5*kilocalories_per_mole/angstrom**2)\n")
        FILE.write("        dummyIndex.append(j)\n")
        FILE.write("        input_positions.append(positions[i])\n")
        FILE.write("\n")
        FILE.write("    integrator = LangevinIntegrator(300*kelvin, 1.0/picoseconds,\n")
        FILE.write("                                      2.0*femtoseconds)\n")
        FILE.write("    context = Context(system, integrator)\n")
        FILE.write("\n")
        FILE.write("    context.setPositions(input_positions)\n")
        FILE.write("\n")
        FILE.write("    run_time = 5.0 # ns\n")
        FILE.write("\n")
        FILE.write("    print('Initial energy:', context.getState(getEnergy=True).getPotentialEnergy())\n")
        FILE.write("    integrator.step(int(run_time * 500000))\n")
        FILE.write("    print('Final energy:', context.getState(getEnergy=True).getPotentialEnergy())\n")
        FILE.write("\n")
        FILE.write("    p = context.getState(getPositions=True).getPositions()[:dummyIndex[0]]\n")
        FILE.write("\n")
        FILE.write("    PDBFile.writeFile(params.topology, p, open('equilibrated.pdb', 'w'))\n")
        FILE.write("\n")
        FILE.write("    return p\n")
        FILE.write("\n")
        FILE.write('"""The actual start of the script"""\n')
        FILE.write("\n")
        FILE.write(f"coords_file = '{structure_file}'\n")
        FILE.write(f"parm_file = '{topology_file}'\n")
        FILE.write('\n')
        #FILE.write("coords = AmberInpcrdFile(coords_file)\n")
        FILE.write("coords = load_file(coords_file)\n")
        #FILE.write("parm = AmberPrmtopFile(parm_file)\n")
        FILE.write("parm = load_file(parm_file)\n")
        FILE.write('\n')
        FILE.write('\n')
        FILE.write("# run the minimisation and equilibration sims\n")
        FILE.write("if not os.path.isfile('equilibrated.pdb'):\n")
        FILE.write("    print('Starting energy minimization')\n")
        FILE.write("    min_pos = run_10k_enmin(parm, coords.positions)\n")
        FILE.write("    print('Done')\n")
        FILE.write("    print('Starting restrained NPT equilibration')\n")
        FILE.write("    npt_pos = run_5ns_NPT_restaint_equil(parm, min_pos)\n")
        FILE.write("    print('Done')\n")
        FILE.write("    print('Starting restrained NVT equilibration')\n")
        FILE.write("    input_pos = run_5ns_NVT_restaint_equil(parm, npt_pos)\n")
        FILE.write("    print('Done')\n")
        FILE.write("else:\n")
        FILE.write("    input_pos = PDBFile('equilibrated.pdb').getPositions()\n")
        FILE.write("\n")
        FILE.write("print('Starting production metadynamics simulation')\n")
        FILE.write("# and now the metadynamics production run\n")
        FILE.write("# prepare system and integrator\n")
        FILE.write("system = parm.createSystem(nonbondedMethod=PME,\n")
        FILE.write("                           nonbondedCutoff=1.0*nanometers,\n")
        FILE.write("                           constraints=HBonds,\n")
        FILE.write("                           rigidWater=True,\n")
        FILE.write("                           ewaldErrorTolerance=0.0005)\n")
        FILE.write("\n")
        FILE.write('\n')
        FILE.write('\n')
        FILE.write(f'lig = [ i for i in range({lig_ids[0]-1}, {lig_ids[-1]})]\n')
        FILE.write(f'p1 = {[ i-1 for i in p1_ids]}\n')
        FILE.write(f'p2 = {[ i-1 for i in p2_ids]}\n')
        FILE.write('\n')
        FILE.write("projection = CustomCentroidBondForce(3, 'distance(g1,g2)*cos(angle(g1,g2,g3))')\n")
        FILE.write('projection.addGroup(lig)\n')
        FILE.write('projection.addGroup(p1)\n')
        FILE.write('projection.addGroup(p2)\n')
        FILE.write('projection.addBond([0,1,2])\n')
        FILE.write('projection.setUsesPeriodicBoundaryConditions(True)\n')
        FILE.write('\n')
        #FILE.write(f'proj = BiasVariable(projection, {(lower_wall-0.2):.2f}*nanometer, {(upper_wall+0.2):.2f}*nanometer, 0.025*nanometer, False, gridWidth = 400)\n')
        FILE.write(f'proj = BiasVariable(projection, {(lower_wall-0.2):.2f}, {(upper_wall+0.2):.2f}, 0.025, False, gridWidth = 200)\n')
        FILE.write('\n')
        FILE.write("extent = CustomCentroidBondForce(3, 'distance(g1,g2)*sin(angle(g1,g2,g3))')\n")
        FILE.write('extent.addGroup(lig)\n')
        FILE.write('extent.addGroup(p1)\n')
        FILE.write('extent.addGroup(p2)\n')
        FILE.write('extent.addBond([0,1,2])\n')
        FILE.write('extent.setUsesPeriodicBoundaryConditions(True)\n')
        FILE.write('\n')
        extent_max = wall_width + wall_buffer + 0.2
#         FILE.write(f'ext = BiasVariable(extent, 0.0*nanometer, {extent_max:.2f}*nanometer, 0.05*nanometer, False, gridWidth = 100)\n')
        FILE.write(f'ext = BiasVariable(extent, 0.0, {extent_max:.2f}, 0.05, False, gridWidth = 200)\n')
        FILE.write('\n')
        FILE.write('# add a flat-bottom restraint\n')
        FILE.write('k1 = 10000*kilojoules_per_mole\n')
        FILE.write('k2 = 1000*kilojoules_per_mole\n')
        FILE.write(f'upper_wall = {upper_wall:.2f}*nanometer\n')
        FILE.write(f'lower_wall = {lower_wall:.2f}\n')
        FILE.write('\n')
        FILE.write('\n')
        FILE.write("upper_wall_rest = CustomCentroidBondForce(3, '(k/2)*max(distance(g1,g2)*cos(angle(g1,g2,g3)) - upper_wall, 0)^2')\n")
        FILE.write('upper_wall_rest.addGroup(lig)\n')
        FILE.write('upper_wall_rest.addGroup(p1)\n')
        FILE.write('upper_wall_rest.addGroup(p2)\n')
        FILE.write('upper_wall_rest.addBond([0,1,2])\n')
        FILE.write("upper_wall_rest.addGlobalParameter('k', k1)\n")
        FILE.write("upper_wall_rest.addGlobalParameter('upper_wall', upper_wall)\n")
        FILE.write('upper_wall_rest.setUsesPeriodicBoundaryConditions(True)\n')
        FILE.write('# dont forget to actually add the force to the system ;)\n')
        FILE.write('system.addForce(upper_wall_rest)\n')
        FILE.write('\n')
        FILE.write('# this is the restraint for the sides of the cone/funnel\n')
        FILE.write('# Ill give a full derivation a bit later in a ipynb\n')
        FILE.write(f'wall_width = {wall_width:.2f}*nanometer\n')
        FILE.write(f'beta_cent = {beta_cent:.2f}\n')
        FILE.write(f's_cent = {s_cent:.2f}*nanometer\n')
        FILE.write(f'wall_buffer = {wall_buffer:.2f}*nanometer\n')
        FILE.write('\n')
        FILE.write("dist_restraint = CustomCentroidBondForce(3, '(k/2)*max(distance(g1,g2)*sin(angle(g1,g2,g3)) - (a/(1+exp(b*(distance(g1,g2)*cos(angle(g1,g2,g3))-c)))+d), 0)^2')\n")
        FILE.write('dist_restraint.addGroup(lig)\n')
        FILE.write('dist_restraint.addGroup(p1)\n')
        FILE.write('dist_restraint.addGroup(p2)\n')
        FILE.write('dist_restraint.addBond([0,1,2])\n')
        FILE.write("dist_restraint.addGlobalParameter('k', k2)\n")
        FILE.write("dist_restraint.addGlobalParameter('a', wall_width)\n")
        FILE.write("dist_restraint.addGlobalParameter('b', beta_cent)\n")
        FILE.write("dist_restraint.addGlobalParameter('c', s_cent)\n")
        FILE.write("dist_restraint.addGlobalParameter('d', wall_buffer)\n")
        FILE.write('dist_restraint.setUsesPeriodicBoundaryConditions(True)\n')
        FILE.write('system.addForce(dist_restraint)\n')
        FILE.write('\n')
        FILE.write('# now the bottom of the funnel\n')
        FILE.write('if lower_wall == 0:\n')
        FILE.write('    #lower_wall = 0.01 didnt work, singularities\n')
        FILE.write('    lower_wall = 0.05 # close enough to zero...\n')
        FILE.write('lower_wall = lower_wall*nanometer\n')
        FILE.write("lower_wall_rest = CustomCentroidBondForce(3, '(k/2)*min(distance(g1,g2)*cos(angle(g1,g2,g3)) - lower_wall, 0)^2')\n")
        FILE.write('lower_wall_rest.addGroup(lig)\n')
        FILE.write('lower_wall_rest.addGroup(p1)\n')
        FILE.write('lower_wall_rest.addGroup(p2)\n')
        FILE.write('lower_wall_rest.addBond([0,1,2])\n')
        FILE.write("lower_wall_rest.addGlobalParameter('k', k1)\n")
        FILE.write("lower_wall_rest.addGlobalParameter('lower_wall', lower_wall)\n")
        FILE.write('lower_wall_rest.setUsesPeriodicBoundaryConditions(True)\n')
        FILE.write('system.addForce(lower_wall_rest)\n')
        FILE.write('\n')
        FILE.write('\n')
        FILE.write('# Set up the simulation.\n')
        FILE.write('\n')
        FILE.write("meta = Metadynamics(system, [proj,ext], 300.0*kelvin, 10.0, 1.5*kilojoules_per_mole, 1000, biasDir = '.', saveFrequency = 250000)\n")
        FILE.write('\n')
        FILE.write('integrator = LangevinIntegrator(300*kelvin, 1.0/picosecond, 0.002*picoseconds)\n')
        FILE.write("platform = Platform.getPlatformByName('CUDA')\n")
        FILE.write("properties = {'CudaPrecision':'mixed'}\n")
        FILE.write('\n')
        FILE.write('simulation = Simulation(parm.topology, system, integrator, platform, properties)\n')
        FILE.write('simulation.context.setPositions(input_pos)\n')
        FILE.write('\n')
        FILE.write('simulation.minimizeEnergy()\n')
        FILE.write('\n')
        #FILE.write("run_time = 1000.0 # ns\n")
        FILE.write(f"run_time = {run_time} # ns\n")
        FILE.write("\n")
        FILE.write("if os.path.isfile('checkpnt.chk'):\n")
        FILE.write("    print('Loading a restart')\n")
        FILE.write("    simulation.loadCheckpoint('checkpnt.chk')\n")
        FILE.write("    shutil.copy('sim_log.csv','old_sim_log.csv')\n")
        FILE.write("    sim_log_file = [ line[:-2] for line in open('sim_log.csv').readlines()]\n")
        FILE.write("    current_steps = int(sim_log_file[-1].split(',')[1])\n")
        FILE.write("    steps = run_time * 500000 - current_steps\n")
        FILE.write("    shutil.copy('COLVAR.npy','old_COLVAR.npy')\n")
        FILE.write("    trj_name = '1_trj.dcd'\n")
        FILE.write("else:\n")
        FILE.write("    steps = run_time * 500000\n")
        FILE.write("    trj_name = 'trj.dcd'\n")
        FILE.write('\n')
        FILE.write("simulation.reporters.append(DCDReporter(trj_name, 50000))\n")
        FILE.write("simulation.reporters.append(CheckpointReporter('checkpnt.chk', 50000))\n")
        FILE.write('simulation.reporters.append(StateDataReporter(\n')
        FILE.write("                           'sim_log.csv', 500000, step=True,\n")
        FILE.write('                           temperature=True,progress=True,\n')
        FILE.write('                           remainingTime=True,speed=True,\n')
        FILE.write("                           totalSteps=steps,separator=','))\n")
        FILE.write('\n')
        FILE.write("def find_bias(filename):\n")
        FILE.write("    return filename.startswith('bias_')\n")
        FILE.write("\n")
        FILE.write("a = np.array([meta.getCollectiveVariables(simulation)])\n")
        FILE.write("for i in range(0,int(steps),1000):\n")
        FILE.write("    if i%10000 == 0:\n")
        FILE.write("        np.save('COLVAR.npy',a)\n")
        FILE.write("    if i%(250000/2) == 0 and i != 0:\n")
        FILE.write("        all_files = os.listdir('.')\n")
        FILE.write("        bias_files = list(filter(find_bias, all_files))\n")
        FILE.write("        latest_bias_file = max(bias_files, key=os.path.getctime)\n")
        FILE.write("        sp.call('cp %s _%s'% (latest_bias_file,latest_bias_file), shell=True)\n")
        FILE.write("    meta.step(simulation, 1000)\n")
        FILE.write("    current_cvs = meta.getCollectiveVariables(simulation)\n")
        FILE.write("    a = np.append(a, [current_cvs], axis=0)\n")
        FILE.write("print('Done')\n")
        FILE.write('\n')
        FILE.write('\n')
