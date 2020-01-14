#Python3

import os
import argparse
import math
import matplotlib.pyplot as plt
import numpy as np
import random as rd


def read_file2(path):
    '''Function that read .db files from a folder'''
    files = []
    #r=root, d=directories, f=files
    for r, d, f in os.walk(path):
        for file in f:
            if '.db' in file:
                files.append(os.path.join(r, file))
    data = []
    for file in files:
        data_file = open(file, 'r', encoding="ISO-8859-1")
        lines = data_file.readlines()
        for i in range(len(lines)):
            if lines[i].find('>') == 0: #Find if the line is the name line
                if lines[i] != lines[0]:
                    data.append(prot)
                prot = []
                l = lines[i].strip() #Remove the \n
                prot.append(l[1:].split()) #Remove the '>' from the name and split
            else:
                prot.append(lines[i].strip().split())
        data_file.close()
    return data


def read_file(db_file):
    '''Function that read a .db file
    Usage: read_file(db_file)
    '''
    data = []
    data_file = open(db_file, 'r', encoding="ISO-8859-1")
    lines = data_file.readlines()
    for i in range(len(lines)):
        if lines[i].find('>') == 0: #Find if the line is the name line
            if lines[i] != lines[0]:
                data.append(prot)
            prot = []
            l = lines[i].strip() #Remove the \n
            prot.append(l[1:].split()) #Remove the '>' from the name and split
        else:
            prot.append(lines[i].strip().split())
    return data


def oneD_seq(data):
    '''Function that takes an array from a .db file and have 2 outputs
    one is a list of Amino acids sequences and other is a list of structure
    sequences, both in 1 dimension
    Usage: oneD_seq(read_file_function_output)
    '''
    aa_seq = []
    str_seq = []
    #Dictionnary of the secondary structure assignement
    dict_struc_numToLet = {'0':'H', '1':'H', '2':'H', '3':'C', '4':'C',
                           '5':'C', '6':'C', '7':'E', '8':'-'}
    for i in range(len(data)):
        amino_seq = ''
        struc_seq = ''
        for j in range(1, len(data[i])):
            amino_seq += data[i][j][0]
            struc_seq += dict_struc_numToLet[data[i][j][2]]
        aa_seq.append(amino_seq)
        str_seq.append(struc_seq)
    return aa_seq, str_seq


def freq_in_seq(seq):
    ''' Function that search for frequence of Amino acids or structures in a list of sequences
    Usage: freq_in_seq(oneD_seq_function_outputs)
    '''
    freq_seq = {}
    tot_length = 0
    for i in range(len(seq)):
        tot_length += len(seq[i])
        for j in range(len(seq[i])):
            #If the observation doesn't exist yet, then we create it in the dictionnary
            if seq[i][j] not in freq_seq:
                freq_seq[seq[i][j]] = 1
            else:
                freq_seq[seq[i][j]] += 1
    for key in freq_seq:
        freq_seq[key] = round(freq_seq[key] * 100 / tot_length, 2)
    return freq_seq


def type_assign(phi1, psi1, phi2, psi2, i3):
    ''' Function that takes phi/psi for residu i+1 and i+2 and amino acid i+2 and return beta turn type
    Usage: type_assign(phi_i+1, psi_i+1, phi_i+2, psi_i+2, i+2_Amino_acid)
    '''
    ecart = 30
    #If amino acid i+2 is a Proline, then special treatment because it is type VI
    if i3 == 'P':
        if (-60 - 45 < phi1 < -60 + 45) and (120 - ecart < psi1 < 120 + ecart) and (-90 - ecart < phi2 < -90 + ecart) and (0 - ecart < psi2 < 0 + ecart):
            return 'typ_6a1'
        elif (-60 - ecart < phi1 < -60 + ecart) and (120 - 45 < psi1 < 120 + 45) and (-90 - ecart < phi2 < -90 + ecart) and (0 - ecart < psi2 < 0 + ecart):
            return 'typ_6a1'
        elif (-60 - ecart < phi1 < -60 + ecart) and (120 - ecart < psi1 < 120 + ecart) and (-90 - 45 < phi2 < -90 + 45) and (0 - ecart < psi2 < 0 + ecart):
            return 'typ_6a1'
        elif (-60 - ecart < phi1 < -60 + ecart) and (120 - ecart < psi1 < 120 + ecart) and (-90 - ecart < phi2 < -90 + ecart) and (0 - 45 < psi2 < 0 + 45):
            return 'typ_6a1'
        elif (-120 - 45 < phi1 < -120 + 45) and (-120 - ecart < psi1 < -120 + ecart) and (-60 - ecart < phi2 < -60 + ecart) and (0 - ecart < psi2 < 0 + ecart):
            return 'typ_6a2'
        elif (-120 - ecart < phi1 < -120 + ecart) and (-120 - 45 < psi1 < -120 + 45) and (-60 - ecart < phi2 < -60 + ecart) and (0 - ecart < psi2 < 0 + ecart):
            return 'typ_6a2'
        elif (-120 - ecart < phi1 < -120 + ecart) and (-120 - ecart < psi1 < -120 + ecart) and (-60 - 45 < phi2 < -60 + 45) and (0 - ecart < psi2 < 0 + ecart):
            return 'typ_6a2'
        elif (-120 - ecart < phi1 < -120 + ecart) and (-120 - ecart < psi1 < -120 + ecart) and (-60 - ecart < phi2 < -60 + ecart) and (0 - 45 < psi2 < 0 + 45):
            return 'typ_6a2'
        elif (-135 - 45 < phi1 < -135 + 45) and (-135 - ecart < psi1 < -135 + ecart) and (-75 - ecart < phi2 < -75 + ecart) and (160 - ecart < psi2 < 160 + ecart):
            return 'typ_6b'
        elif (-135 - ecart < phi1 < -135 + ecart) and (-135 - 45 < psi1 < -135 + 45) and (-75 - ecart < phi2 < -75 + ecart) and (160 - ecart < psi2 < 160 + ecart):
            return 'typ_6b'
        elif (-135 - ecart < phi1 < -135 + ecart) and (-135 - ecart < psi1 < -135 + ecart) and (-75 - 45 < phi2 < -75 + 45) and (160 - ecart < psi2 < 160 + ecart):
            return 'typ_6b'
        elif (-135 - ecart < phi1 < -135 + ecart) and (-135 - ecart < psi1 < -135 + ecart) and (-75 - ecart < phi2 < -75 + ecart) and (160 - 45 < psi2 < 160 + 45):
            return 'typ_6b'
    elif (-60 - 45 < phi1 < -60 + 45) and (120 - ecart < psi1 < 120 + ecart) and (-90 - ecart < phi2 < -90 + ecart) and (0 - ecart < psi2 < 0 + ecart):
        return 'typ_1'
    elif (-60 - ecart < phi1 < -60 + ecart) and (120 - 45 < psi1 < 120 + 45) and (-90 - ecart < phi2 < -90 + ecart) and (0 - ecart < psi2 < 0 + ecart):
        return 'typ_1'
    elif (-60 - ecart < phi1 < -60 + ecart) and (120 - ecart < psi1 < 120 + ecart) and (-90 - 45 < phi2 < -90 + 45) and (0 - ecart < psi2 < 0 + ecart):
        return 'typ_1'
    elif (-60 - ecart < phi1 < -60 + ecart) and (120 - ecart < psi1 < 120 + ecart) and (-90 - ecart < phi2 < -90 + ecart) and (0 - 45 < psi2 < 0 + 45):
        return 'typ_1'
    elif (60 - 45 < phi1 < 60 + 45) and (-120 - ecart < psi1 < -120 + ecart) and (90 - ecart < phi2 < 90 + ecart) and (0 - ecart < psi2 < 0 + ecart):
        return 'typ_1p'
    elif (60 - ecart < phi1 < 60 + ecart) and (-120 - 45 < psi1 < -120 + 45) and (90 - ecart < phi2 < 90 + ecart) and (0 - ecart < psi2 < 0 + ecart):
        return 'typ_1p'
    elif (60 - ecart < phi1 < 60 + ecart) and (-120 - ecart < psi1 < -120 + ecart) and (90 - 45 < phi2 < 90 + 45) and (0 - ecart < psi2 < 0 + ecart):
        return 'typ_1p'
    elif (60 - ecart < phi1 < 60 + ecart) and (-120 - ecart < psi1 < -120 + ecart) and (90 - ecart < phi2 < 90 + ecart) and (0 - 45 < psi2 < 0 + 45):
        return 'typ_1p'
    elif (-60 - 45 < phi1 < -60 + 45) and (120 - ecart < psi1 < 120 + ecart) and (80 - ecart < phi2 < 80 + ecart) and (0 - ecart < psi2 < 0 + ecart):
        return 'typ_2'
    elif (-60 - ecart < phi1 < -60 + ecart) and (120 - 45 < psi1 < 120 + 45) and (80 - ecart < phi2 < 80 + ecart) and (0 - ecart < psi2 < 0 + ecart):
        return 'typ_2'
    elif (-60 - ecart < phi1 < -60 + ecart) and (120 - ecart < psi1 < 120 + ecart) and (80 - 45 < phi2 < 80 + 45) and (0 - ecart < psi2 < 0 + ecart):
        return 'typ_2'
    elif (-60 - ecart < phi1 < -60 + ecart) and (120 - ecart < psi1 < 120 + ecart) and (80 - ecart < phi2 < 80 + ecart) and (0 - 45 < psi2 < 0 + 45):
        return 'typ_2'
    elif (60 - 45 < phi1 < 60 + 45) and (-120 - ecart < psi1 < -120 + ecart) and (-80 - ecart < phi2 < -80 + ecart) and (0 - ecart < psi2 < 0 + ecart):
        return 'typ_2p'
    elif (60 - ecart < phi1 < 60 + ecart) and (-120 - 45 < psi1 < -120 + 45) and (-80 - ecart < phi2 < -80 + ecart) and (0 - ecart < psi2 < 0 + ecart):
        return 'typ_2p'
    elif (60 - ecart < phi1 < 60 + ecart) and (-120 - ecart < psi1 < -120 + ecart) and (-80 - 45 < phi2 < -80 + 45) and (0 - ecart < psi2 < 0 + ecart):
        return 'typ_2p'
    elif (60 - ecart < phi1 < 60 + ecart) and (-120 - ecart < psi1 < -120 + ecart) and (-80 - ecart < phi2 < -80 + ecart) and (0 - 45 < psi2 < 0 + 45):
        return 'typ_2p'
    elif (-60 - 45 < phi1 < -60 + 45) and (-30 - ecart < psi1 < -30 + ecart) and (-120 - ecart < phi2 < -120 + ecart) and (120 - ecart < psi2 < 120 + ecart):
        return 'typ_8'
    elif (-60 - ecart < phi1 < -60 + ecart) and (-30 - 45 < psi1 < -30 + 45) and (-120 - ecart < phi2 < -120 + ecart) and (120 - ecart < psi2 < 120 + ecart):
        return 'typ_8'
    elif (-60 - ecart < phi1 < -60 + ecart) and (-30 - ecart < psi1 < -30 + ecart) and (-120 - 45 < phi2 < -120 + 45) and (120 - ecart < psi2 < 120 + ecart):
        return 'typ_8'
    elif (-60 - ecart < phi1 < -60 + ecart) and (-30 - ecart < psi1 < -30 + ecart) and (-120 - ecart < phi2 < -120 + ecart) and (120 - 45 < psi2 < 120 + 45):
        return 'typ_8'
    elif (-120 - 45 < phi1 < -120 + 45) and (130 - ecart < psi1 < 130 + ecart) and (55 - ecart < phi2 < 55 + ecart) and (41 - ecart < psi2 < 41 + ecart):
        return 'typ_4-1'
    elif (-120 - ecart < phi1 < -120 + ecart) and (130 - 45 < psi1 < 130 + 45) and (55 - ecart < phi2 < 55 + ecart) and (41 - ecart < psi2 < 41 + ecart):
        return 'typ_4-1'
    elif (-120 - ecart < phi1 < -120 + ecart) and (130 - ecart < psi1 < 130 + ecart) and (55 - 45 < phi2 < 55 + 45) and (41 - ecart < psi2 < 41 + ecart):
        return 'typ_4-1'
    elif (-120 - ecart < phi1 < -120 + ecart) and (130 - ecart < psi1 < 130 + ecart) and (55 - ecart < phi2 < 55 + ecart) and (41 - 45 < psi2 < 41 + 45):
        return 'typ_4-1'
    elif (-85 - 45 < phi1 < -85 + 45) and (-15 - ecart < psi1 < -15 + ecart) and (-125 - ecart < phi2 < -125 + ecart) and (55 - ecart < psi2 < 55 + ecart):
        return 'typ_4-2'
    elif (-85 - ecart < phi1 < -85 + ecart) and (-15 - 45 < psi1 < -15 + 45) and (-125 - ecart < phi2 < -125 + ecart) and (55 - ecart < psi2 < 55 + ecart):
        return 'typ_4-2'
    elif (-85 - ecart < phi1 < -85 + ecart) and (-15 - ecart < psi1 < -15 + ecart) and (-125 - 45 < phi2 < -125 + 45) and (55 - ecart < psi2 < 55 + ecart):
        return 'typ_4-2'
    elif (-85 - ecart < phi1 < -85 + ecart) and (-15 - ecart < psi1 < -15 + ecart) and (-125 - ecart < phi2 < -125 + ecart) and (55 - 45 < psi2 < 55 + 45):
        return 'typ_4-2'
    elif (-71 - 45 < phi1 < -71 + 45) and (-30 - ecart < psi1 < -30 + ecart) and (-72 - ecart < phi2 < -72 + ecart) and (-47 - ecart < psi2 < -47 + ecart):
        return 'typ_4-3'
    elif (-71 - ecart < phi1 < -71 + ecart) and (-30 - 45 < psi1 < -30 + 45) and (-72 - ecart < phi2 < -72 + ecart) and (-47 - ecart < psi2 < -47 + ecart):
        return 'typ_4-3'
    elif (-71 - ecart < phi1 < -71 + ecart) and (-30 - ecart < psi1 < -30 + ecart) and (-72 - 45 < phi2 < -72 + 45) and (-47 - ecart < psi2 < -47 + ecart):
        return 'typ_4-3'
    elif (-71 - ecart < phi1 < -71 + ecart) and (-30 - ecart < psi1 < -30 + ecart) and (-72 - ecart < phi2 < -72 + ecart) and (-47 - 45 < psi2 < -47 + 45):
        return 'typ_4-3'
    elif (-97 - 45 < phi1 < -97 + 45) and (-2 - ecart < psi1 < -2 + ecart) and (-117 - ecart < phi2 < -117 + ecart) and (-11 - ecart < psi2 < -11 + ecart):
        return 'typ_4-4'
    elif (-97 - ecart < phi1 < -97 + ecart) and (-2 - 45 < psi1 < -2 + 45) and (-117 - ecart < phi2 < -117 + ecart) and (-11 - ecart < psi2 < -11 + ecart):
        return 'typ_4-4'
    elif (-97 - ecart < phi1 < -97 + ecart) and (-2 - ecart < psi1 < -2 + ecart) and (-117 - 45 < phi2 < -117 + 45) and (-11 - ecart < psi2 < -11 + ecart):
        return 'typ_4-4'
    elif (-97 - ecart < phi1 < -97 + ecart) and (-2 - ecart < psi1 < -2 + ecart) and (-117 - ecart < phi2 < -117 + ecart) and (-11 - 45 < psi2 < -11 + 45):
        return 'typ_4-4'
    else:
        return 'typ_4misc'


def kohonen_files(struc_seq, data):
    '''Return a list with phi.spi angles of the residues and a list with the turn type of the turns
    Usage: kohonen_files(oneD_seq_Function_output_struc_seq, read_file_function_output)
    '''
    type_file = []
    angles_file = []
    for i in range(len(struc_seq)):
        for j in range(len(struc_seq[i])-3):
            #Distance between residue i and i+3
            dist = math.sqrt((float(data[i][j+4][6]) - float(data[i][j+1][6]))**2 + (float(data[i][j+4][7]) - float(data[i][j+1][7]))**2 + (float(data[i][j+4][8]) - float(data[i][j+1][8]))**2)
            #Condition for 4 residues to be a beta turn
            if (dist <= 7) and (str(struc_seq[i][j:j+4]) != 'EEEE') and str(struc_seq[i][j+1] + struc_seq[i][j+2]) != 'HH':
                phi1 = float(data[i][j+2][3])
                psi1 = float(data[i][j+2][4])
                phi2 = float(data[i][j+3][3])
                psi2 = float(data[i][j+3][4])
                angles_file.append([phi1, psi1, phi2, psi2])
                type_file.append(type_assign(phi1, psi1, phi2, psi2, data[i][j+3][0]))
    return type_file, angles_file


def count_type(typ_list):
    ''' Function that takes a list of type and return the frequencies of each.
    Usage: count_type(kohonen_files_function_output_type_file)
    '''
    numb_type = {}
    length = len(typ_list)
    #Cout number of types as frequences (percent)
    numb_type['typ_1'] = round(typ_list.count('typ_1') * 100 / length, 2)
    numb_type['typ_1p'] = round(typ_list.count('typ_1p') * 100 / length, 2)
    numb_type['typ_2'] = round(typ_list.count('typ_2') * 100 / length, 2)
    numb_type['typ_2p'] = round(typ_list.count('typ_2p') * 100 / length, 2)
    numb_type['typ_4-1'] = round(typ_list.count('typ_4-1') * 100 / length, 2)
    numb_type['typ_4-2'] = round(typ_list.count('typ_4-2') * 100 / length, 2)
    numb_type['typ_4-3'] = round(typ_list.count('typ_4-3') * 100 / length, 2)
    numb_type['typ_4-4'] = round(typ_list.count('typ_4-4') * 100 / length, 2)
    numb_type['typ_4misc'] = round(typ_list.count('typ_4misc') * 100 / length, 2)
    numb_type['typ_6a1'] = round(typ_list.count('typ_6a1') * 100 / length, 2)
    numb_type['typ_6a2'] = round(typ_list.count('typ_6a2') * 100 / length, 2)
    numb_type['typ_6b'] = round(typ_list.count('typ_6b') * 100 / length, 2)
    numb_type['typ_8'] = round(typ_list.count('typ_8') * 100 / length, 2)
    return numb_type


def distance(val, val_node):
    ''' Function that takes the values of a node and the values of angles and return the distance between them
    Usage: distance(angles_values, node_values)
    '''
    dist = abs(val_node - val)
    if dist > 180:
        dist = 360 - dist
    return dist


def plot_map(kohonen_learn, types, title, save):
    """Function to plot ramachandran from a Koonen map
    Title is the title we want on the plot
    Save is the name of the file that will be saved
    Usage: plot_map(kohonen_map, title_of_the_plot, name_of_the_seved_plot)
    """
    plt.rcParams["figure.figsize"] = [17,10]
    fig, axs = plt.subplots(4,4, sharex = True, sharey = True)
    lign = 0
    col = 0
    n = 0
    for i in range(len(kohonen_learn)):
        for j in range(len(kohonen_learn[i])):
            axs[lign, col].plot(kohonen_learn[i][j][0], kohonen_learn[i][j][1], 'ro', label = 'Residu i+1')
            axs[lign, col].plot(kohonen_learn[i][j][2], kohonen_learn[i][j][3], 'bo', label = 'Residu i+2')
            axs[lign, col].plot([-180,180], [0,0], 'k')
            axs[lign, col].plot([0,0], [-180,180], 'k')
            axs[lign, col].set_xlim(-180, 180)
            axs[lign, col].set_ylim(-180, 180)
            axs[lign, col].set_title(types[n])
            axs[lign, col].set_aspect('equal',adjustable='box')
            col += 1
            n += 1
        col = 0
        lign += 1
    for ax in axs.flat:
        ax.set(xlabel='Phi', ylabel='Psi')
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center')
    fig.suptitle(title)
    plt.rcParams["figure.figsize"] = [16,9]
    plt.savefig('../results/shcullpdb_pc20_res2.5_R1.0_d050827_chains2722/' + save)
    plt.show()
    return


def kohonen(size, angle_files, nb_loop):
    '''Kohonen maps
    size = (X,Y), with X and Y as integers
    angles_files is the lsit containing the angles phi/psi for residues i+1 and i+2
    nb_loop is the number of learning iteration
    Usage: kohonen(size_of_the_neural_map, angle_list, iteration_number)
    '''
    #Shuffleing
    rd.shuffle(angle_files)

    nb_neuron = size[0] * size[1]
    size.append(4) #Number of angles : 2 for residue i+1 and 2 for i+2
    koho_map = np.zeros(size, dtype=float)
    #Initialisation
    for i in range(len(koho_map)):
        for j in range(len(koho_map[i])):
            val_init = rd.choice(angle_files)
            koho_map[i][j] = val_init

    #Analyse type of the neurons of the map
    type_carte = []
    for i in range(len(koho_map)):
        for j in range(len(koho_map[i])):
            type_carte.append(type_assign(koho_map[i][j][0], koho_map[i][j][1], koho_map[i][j][2], koho_map[i][j][3], 'A')) #We put amino acid i+2 to A here to not make the function has a error
    #Generate initial map
    plot_map(koho_map, type_carte, 'Carte pre apprentissage', 'carte_Pre_Apprentissage.png')

    a0 = 0.8 #Alpha init
    n0 = 4 #Beta init
    for loop in range(nb_loop):
        val = rd.choice(angle_files)
        #Search Winner using distance
        d = 200 #Initialise a large distance
        for i in range(len(koho_map)):
            for j in range(len(koho_map[i])):
                #Calculate the total distance
                dist = math.sqrt(distance(val[0], koho_map[i][j][0])**2 + distance(val[1], koho_map[i][j][1])**2 + distance(val[2], koho_map[i][j][2])**2 + distance(val[3], koho_map[i][j][3])**2)
                #If the distance is shorter, then new winner with new best distance
                if dist < d:
                    d = dist
                    winner = [i,j]

        #Find first neighbors nodes (8 nodes)
        neighbors = [[winner[0]-1, winner[1]-1], [winner[0]-1, winner[1]], [winner[0]-1, winner[1]+1], [winner[0], winner[1]-1], [winner[0], winner[1]+1], [winner[0]+1, winner[1]-1], [winner[0]+1, winner[1]], [winner[0]+1, winner[1]+1]]
        #Create distance array between winner and neighbors
        dist_node_to_Wnode = []
        for i in range(len(neighbors)):
            dist_node_to_Wnode.append(math.sqrt((winner[0] - neighbors[i][0])**2 + (winner[1] - neighbors[i][1])**2))
        #Loop to close the map and allow nodes to be linked
        for i in range(len(neighbors)):
            if neighbors[i][0] > len(koho_map)-1:
                neighbors[i][0] = 0
            if neighbors[i][1] > len(koho_map[0])-1:
                neighbors[i][1] = 0
            if neighbors[i][0] == -1:
                neighbors[i][0] = len(koho_map)-1
            if neighbors[i][1] == -1:
                neighbors[i][1] = len(koho_map[0])-1

        coeff_a = (a0) / (1 + (loop / nb_loop))
        nu = (n0) / (1 + (loop / nb_loop))
        coeff_b = np.exp(- (0)**2 / (2*(nu**2))) #Distance between winner and winner is 0, then the coeff_b for the winner is 1
        gamma = coeff_a * coeff_b
        #Modification of the value of the winner
        koho_map[winner[0],winner[1]] = koho_map[winner[0], winner[1]] + (val - koho_map[winner[0],winner[1]]) * gamma

        #Change nodes values (diffusion) in function of tdistance to the winner
        for i in range(len(neighbors)):
            coeff_b = np.exp(- (dist_node_to_Wnode[i])**2 / (2*(nu**2)))
            gamma = coeff_a * coeff_b
            koho_map[neighbors[i][0],neighbors[i][1]] = koho_map[neighbors[i][0],neighbors[i][1]] + (val - koho_map[neighbors[i][0],neighbors[i][1]]) * gamma
    return koho_map


def write_dict_on_file(freq_AA, file_name, header):
    ''' Function that takes a list and write it in a file
    Usage: write_dict_on_file(list_to_write, name_of_the_file, header_of_the_file)
    '''
    if not os.path.exists('../results/shcullpdb_pc20_res2.5_R1.0_d050827_chains2722'):
        os.makedirs('../results/shcullpdb_pc20_res2.5_R1.0_d050827_chains2722')
    f = open('../results/shcullpdb_pc20_res2.5_R1.0_d050827_chains2722/' + file_name, 'w')
    f.write(header + '\n')
    for key, value in freq_AA:
        f.write(str(key) + '\t' + str(value) + '\n')
    return


def main():
    """Main function"""
    #Arguments to run the program
    parser = argparse.ArgumentParser(description='Kohonen map & Beta Turn')
    parser.add_argument('-f', metavar='File', type=str,help='db File', required=True)
    args = parser.parse_args()
    #Import arguments as variables
    db_file = args.f
    data = read_file(args.f)

    #Dictionnary to transform sequences in Amino Acids (AA) and strucutres (HEC)
    data_prot_AA, data_prot_Struc = oneD_seq(data)

    #Frequency of AA
    freq_AA_seq = freq_in_seq(data_prot_AA)
    #Frequency of HEC
    freq_Struc_seq = freq_in_seq(data_prot_Struc)
    #2 lists: one with the angles and one wiht the types
    typ_file, angle_file = kohonen_files(data_prot_Struc, data)
    #Count types
    numb_type = count_type(typ_file)

    #Write on files the different frequencies
    write_dict_on_file(freq_AA_seq.items(), 'aa_freq.txt', 'Amino_Acids' + '\t' + 'Frequency_(%)')
    write_dict_on_file(freq_Struc_seq.items(), 'struct_freq.txt', 'Secondary_Structure' + '\t' + 'Frequency_(%)')
    write_dict_on_file(numb_type.items(), 'type_freq.txt', 'Type' + '\t' + 'Frequency_(%)')

    #Initialise the kohonen map and learning
    kohonen_learn = kohonen([4,4], angle_file, 200)

    #Analyse type of the neurons of the map
    type_carte = []
    for i in range(len(kohonen_learn)):
        for j in range(len(kohonen_learn[i])):
            type_carte.append(type_assign(kohonen_learn[i][j][0], kohonen_learn[i][j][1], kohonen_learn[i][j][2], kohonen_learn[i][j][3], 'A')) #We put amino acid i+2 to A here to not make the function has a error

    #Plot kohonen map
    plot_map(kohonen_learn, type_carte, 'Carte post apprentissage', 'carte_Post_Apprentissage.png')

if __name__ == '__main__':
    main()
