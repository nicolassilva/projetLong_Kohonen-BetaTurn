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
    '''Function that read a .db file'''
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
    sequences, both in 1 dimension'''
    aa_seq = []
    str_seq = []
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
    ''' Function that search for frequence of Amino acids or structures or coil
    in a list of sequences'''
    freq_seq = {}
    tot_length = 0
    for i in range(len(seq)):
        tot_length += len(seq[i])
        for j in range(len(seq[i])):
            if seq[i][j] not in freq_seq:
                freq_seq[seq[i][j]] = 1
            else:
                freq_seq[seq[i][j]] += 1
    for key in freq_seq:
        freq_seq[key] = round(freq_seq[key] * 100 / tot_length, 2)
    return freq_seq


def kohonen_files(struc_seq, data, ecart):
    '''Return a list with phi.spi angles of the residues and a list with the turn type of the turns
    ecart is the error marge between the value of the angle and the acceptable value of the turn type'''
    type_file = []
    angles_file = []
    for i in range(len(struc_seq)):
        for j in range(len(struc_seq[i])-3):
            #Distance between residue i and i+3
            dist = math.sqrt((float(data[i][j+4][6]) - float(data[i][j+1][6]))**2 + (float(data[i][j+4][7]) - float(data[i][j+1][7]))**2 + (float(data[i][j+4][8]) - float(data[i][j+1][8]))**2)
            #print(dist)
            if struc_seq[i][j] != struc_seq[i][j+1] or struc_seq[i][j] != struc_seq[i][j+2] or struc_seq[i][j] != struc_seq[i][j+3] and struc_seq[i][j+1] != 'H' and struc_seq[i][j+2] != 'H' and dist < 7 and struc_seq[i][j] != struc_seq[i][j+1] != struc_seq[i][j+2] != struc_seq[i][j+3] != 'E':
                phi1 = float(data[i][j+2][3])
                psi1 = float(data[i][j+2][4])
                phi2 = float(data[i][j+3][3])
                psi2 = float(data[i][j+3][4])
                angles_file.append([phi1, psi1, phi2, psi2])
                if (-60 - ecart < phi1 < -60 + ecart) and (-30 - ecart < psi1 < -30 + ecart) and (-90 - ecart < phi2 < -90 + ecart) and (0 - ecart < psi2 < 0 + ecart):
                    type_file.append('typ_1')
                elif (60 - ecart < phi1 < 60 + ecart) and (30 - ecart < psi1 < 30 + ecart) and (90 - ecart < phi2 < 90 + ecart) and (0 - ecart < psi2 < 0 + ecart):
                    type_file.append('typ_1p')
                elif (-60 - ecart < phi1 < -60 + ecart) and (120 - ecart < psi1 < 120 + ecart) and (80 - ecart < phi2 < 80 + ecart) and (0 - ecart < psi2 < 0 + ecart):
                    type_file.append('typ_2')
                elif (60 - ecart < phi1 < 60 + ecart) and (-120 - ecart < psi1 < -120 + ecart) and (-80 - ecart < phi2 < -80 + ecart) and (0 - ecart < psi2 < 0 + ecart):
                    type_file.append('typ_2p')
                elif (-60 - ecart < phi1 < -60 + ecart) and (120 - ecart < psi1 < 120 + ecart) and (-90 - ecart < phi2 < -90 + ecart) and (0 - ecart < psi2 < 0 + ecart):
                    type_file.append('typ_6a1')
                elif (-120 - ecart < phi1 < -120 + ecart) and (-120 - ecart < psi1 < -120 + ecart) and (-60 - ecart < phi2 < -60 + ecart) and (0 - ecart < psi2 < 0 + ecart):
                    type_file.append('typ_6a2')
                elif (-135 - ecart < phi1 < -135 + ecart) and (-135 - ecart < psi1 < -135 + ecart) and (-75 - ecart < phi2 < -75 + ecart) and (160 - ecart < psi2 < 160 + ecart):
                    type_file.append('typ_6b')
                elif (-120 - ecart < phi1 < -120 + ecart) and (130 - ecart < psi1 < 130 + ecart) and (55 - ecart < phi2 < 55 + ecart) and (41 - ecart < psi2 < 41 + ecart):
                    type_file.append('typ_4-1')
                elif (-85 - ecart < phi1 < -85 + ecart) and (-15 - ecart < psi1 < -15 + ecart) and (-125 - ecart < phi2 < -125 + ecart) and (55 - ecart < psi2 < 55 + ecart):
                    type_file.append('typ_4-2')
                elif (-71 - ecart < phi1 < -71 + ecart) and (-30 - ecart < psi1 < -30 + ecart) and (-72 - ecart < phi2 < -72 + ecart) and (-47 - ecart < psi2 < -47 + ecart):
                    type_file.append('typ_4-3')
                elif (-97 - ecart < phi1 < -97 + ecart) and (-2 - ecart < psi1 < -2 + ecart) and (-117 - ecart < phi2 < -117 + ecart) and (-11 - ecart < psi2 < -11 + ecart):
                    type_file.append('typ_4-4')
                elif (-60 - ecart < phi1 < -60 + ecart) and (-30 - ecart < psi1 < -30 + ecart) and (-120 - ecart < phi2 < -120 + ecart) and (120 - ecart < psi2 < 120 + ecart):
                    type_file.append('typ_8')
                else:
                    type_file.append('typ_4misc')
    return type_file, angles_file


def count_type(typ_list):
    numb_type = {}
    length = len(typ_list)
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


def ramachandran(angle_file, typ_file):
    """Function that takes angles phi/psi for residue i+1 and i+2 and turn type and plot ramachandran graph of the data"""
    dict_col = {'typ_1':'gold', 'typ_1p':'darkkhaki', 'typ_2':'darkblue', 'typ_2p':'violet', 'typ_4-1':'green', 'typ_4-2':'darkgreen', 'typ_4-3':'limegreen', 'typ_4-4':'greenyellow', 'typ_6a1':'darkgrey', 'typ_6a2':'brown', 'typ_6b':'lightsalmon', 'typ_8':'red'}
    col = ['gold', 'darkkhaki', 'darkblue', 'violet',  'green', 'darkgreen', 'limegreen', 'greenyellow', 'darkgrey', 'brown', 'lightsalmon', 'red']
    typ_name = ['typ_1', 'typ_1p', 'typ_2', 'typ_2p', 'typ_4-1', 'typ_4-2', 'typ_4-3', 'typ_4-4', 'typ_6a1', 'typ_6a2', 'typ_6b', 'typ_8']

    size_axes = 180
    plt.figure(1)
    plt.subplot(121)
    plt.title('Ramachandran plot i+1')
    plt.xlabel("Phi")
    plt.ylabel("Psi")
    plt.plot([-size_axes,size_axes], [0,0], 'k')
    plt.plot([0,0], [-size_axes,size_axes], 'k')
    plt.gca().set_aspect('equal',adjustable='box')
    plt.xlim((-size_axes, size_axes))
    plt.ylim((-size_axes, size_axes))
    plt.subplot(122)
    plt.title('Ramachandran plot i+2')
    plt.xlabel("Phi")
    plt.ylabel("Psi")
    plt.plot([-size_axes,size_axes], [0,0], 'k')
    plt.plot([0,0], [-size_axes,size_axes], 'k')
    plt.gca().set_aspect('equal',adjustable='box')
    plt.xlim((-size_axes, size_axes))
    plt.ylim((-size_axes, size_axes))

    for i in range(len(angle_file)):
        plt.subplot(121)
        plt.plot(angle_file[i][0], angle_file[i][1], 'o', c=dict_col[str(typ_file[i])], label=str(typ_file[i]))
        plt.subplot(122)
        plt.plot(angle_file[i][2], angle_file[i][3], 'o', c=dict_col[str(typ_file[i])], label=str(typ_file[i]))

    plt.text(-680, -280, "Type 1", c='gold')
    plt.text(-680, -320, "Type 1'", c='darkkhaki')
    plt.text(-530, -280, "Type 2", c='darkblue')
    plt.text(-530, -320, "Type 2'", c='violet')
    plt.text(-380, -280, "Type 5a1", c='green')
    plt.text(-380, -320, "Type 5a2", c='darkgreen')
    plt.text(-230, -280, "Type 5b", c='limegreen')
    plt.text(-230, -320, "Type 4_1", c='greenyellow')
    plt.text(-80, -280, "Type 4_2", c='darkgrey')
    plt.text(-80, -320, "Type 4_3", c='brown')
    plt.text(70, -280, "Type 4_8", c='lightsalmon')
    plt.text(70, -320, "Type 8", c='red')
    plt.show()
    return


def distance(val, val_node):
    dist = abs(val_node - val)
    if dist > 180:
        dist = 360 - dist
    return dist


def kohonen(size, angle_files, nb_loop):
    '''Kohonen maps
    size = (X,Y), with X and Y as integers
    col_phi and col_psi are the column in the angles files: col 0 and 1 are the phi/psi for residue i+1 and col 2 and 3 are the phi/psi for residue i+2
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
    print(koho_map)
    #Search Winner using distance
    a0 = 0.8 #Alpha init
    n0 = 0.5 #Beta init
    for loop in range(1):
        val = rd.choice(angle_files)
        #print(val)
        d = 200 #Initialise a large distance
        #
        for i in range(len(koho_map)):
            for j in range(len(koho_map[i])):
                #Calculate the total distance
                dist = math.sqrt(distance(val[0], koho_map[i][j][0])**2 + distance(val[1], koho_map[i][j][1])**2 + distance(val[2], koho_map[i][j][2])**2 + distance(val[3], koho_map[i][j][3])**2)
                #If the distance is shorter, then new winner with new best distance
                if dist < d:
                    d = dist
                    winner = [i,j]
        print(winner)

        #Find first neighbors nodes (8 nodes)
        neighbors = [[winner[0]-1, winner[1]-1], [winner[0]-1, winner[1]], [winner[0]-1, winner[1]+1], [winner[0], winner[1]-1], [winner[0], winner[1]+1], [winner[0]+1, winner[1]-1], [winner[0]+1, winner[1]], [winner[0]+1, winner[1]+1]]
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
        #print(neighbors)
        #print(koho_map)

        coeff_a = (a0) / (1 + (loop / nb_loop))
        coeff_b = 4
        gamma = coeff_a * coeff_b
        #Modification of the value of the winner
        koho_map[winner[0],winner[1]] = koho_map[winner[0],winner[1]] + (val - koho_map[winner[0],winner[1]]) * gamma
        #print(koho_map)

        #Change nodes values (diffusion)
        nu = (n0) / (1 + (loop / nb_loop))
        for i in range(len(neighbors)):
            dist_node_to_Wnode = math.sqrt((winner[0] - neighbors[i][0])**2 + (winner[1] - neighbors[i][1])**2)
            coeff_b = np.exp(- (dist_node_to_Wnode)**2 / (2*(nu**2)))
            #print(coeff_b)
            gamma = coeff_a * coeff_b
            #print(koho_map[neighbors[i][0],neighbors[i][1]])
            koho_map[neighbors[i][0],neighbors[i][1]] = koho_map[neighbors[i][0],neighbors[i][1]] + (val - koho_map[neighbors[i][0],neighbors[i][1]]) * gamma
            #print(gamma)
            #print(koho_map[neighbors[i][0],neighbors[i][1]])
    return koho_map


def write_dict_on_file(freq_AA, file_name, header):
    if not os.path.exists('../results/shcullpdb_pc20_res2.5_R1.0_d050827_chains2722'):
        os.makedirs('../results/shcullpdb_pc20_res2.5_R1.0_d050827_chains2722')
    f = open('../results/shcullpdb_pc20_res2.5_R1.0_d050827_chains2722/' + file_name, 'w')
    f.write(header + '\n')
    for key, value in freq_AA:
        f.write(str(key) + '\t' + str(value) + '\n')
    return


def plot_map(kohonen_learn, X, title):
    """Function to plot ramachandran from a Koonen map
    X is the subplot number
    """
    plt.subplot(X)
    x = []
    y = []
    for i in range(len(kohonen_learn)):
        for j in range(len(kohonen_learn[i])):
            x.append(kohonen_learn[i][j][0])
            y.append(kohonen_learn[i][j][1])
    plt.plot(x, y, 'o')
    plt.plot([-180,180], [0,0], 'k')
    plt.plot([0,0], [-180,180], 'k')
    plt.xlim(-180, 180)
    plt.xlabel('Phi')
    plt.ylim(-180, 180)
    plt.ylabel('Psi')
    plt.title(title)
    plt.gca().set_aspect('equal',adjustable='box')
    #plt.show()
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
    typ_file, angle_file = kohonen_files(data_prot_Struc, data, 30)
    #Count types
    numb_type = count_type(typ_file)

    #Write on files the different frequencies
    write_dict_on_file(freq_AA_seq.items(), 'aa_freq.txt', 'Amino_Acids' + '\t' + 'Frequency_(%)')
    write_dict_on_file(freq_Struc_seq.items(), 'struct_freq.txt', 'Secondary_Structure' + '\t' + 'Frequency_(%)')
    write_dict_on_file(numb_type.items(), 'type_freq.txt', 'Type' + '\t' + 'Frequency_(%)')

    #Ramachandran plot of the data
    #ramachandran(angle_file, typ_file)

    #Initialise the kohonen map and learning
    kohonen_learn = kohonen([4,4], angle_file, 25)
    print(kohonen_learn)

    '''#Plot kohonen map
    plt.figure(1)
    plot_map(kohonen_learn_2, 122, 'Résidues i+2')
    plt.show()'''

if __name__ == '__main__':
    main()
