import cleanup
import classify
import wavyness
import somaVolume
import kinkPositions
import beads
import os
import time
import csv
from tqdm import tqdm
from winsound import MessageBeep
import utility
#np.set_printoptions(precision=2, suppress=True)

start = time.time()

#Parameters
cleanup_tree = True
scale = (0.223, 0.223, 0.3) #um
length_threshold = 3 #um
angle_threshold = 145 #degree
window_size_linear_regression = 6.0 #um
window_size_maximum_supression = 5.0 #um
plm_outgrowth_suppression = 5.0 #um length of plm mainbranch at the end that is suppressed for outgrowth as the neuron tracing algorithm often makes mistakes at the end
neurontypes = ['ALM', 'PLM'] #list containing 'ALM', 'PLM' or both
soma_volume_singleseed = True
visualize_somavolumes = False
root = 'E:/NeuronMorphologyQuantificationData/'
infolder = 'trees/'
infolder_tif = 'images/'


outfolder_cleanedtrees = 'cleanedtrees/'
outfolder_classifiedtrees = 'classifiedtrees/'
outfolder_wavytrees = 'wavytrees/'
outfolder_gifs = 'somavolume_gifs/'
now = time.strftime("%Y-%m-%d_%H-%M-%S")
measurement_filename = now + '_IndividualMeasurements.csv'
measurement_filename2 = now + '_SummaryMeasurements.csv'
parameters_filename = now + '_Parameters.txt'

#Write the parameter file
with open(root + parameters_filename, 'w') as f:
    lj = 30
    f.write('Clean tree:'.ljust(lj) + '{}\n'.format(cleanup_tree))
    f.write('Scale:'.ljust(lj) + '{}\n'.format(scale))
    f.write('Length threshold:'.ljust(lj) + '{}\n'.format(length_threshold))
    f.write('Angle Threshold:'.ljust(lj) + '{}\n'.format(angle_threshold))
    f.write('Window lin reg:'.ljust(lj) + '{}\n'.format(window_size_linear_regression))    
    f.write('Window max supp:'.ljust(lj) + '{}\n'.format(window_size_maximum_supression))
    f.write('PLM outgrowth supp:'.ljust(lj) + '{}\n'.format(plm_outgrowth_suppression))
    f.write('Soma volume single seed:'.ljust(lj) + '{}\n'.format(soma_volume_singleseed))
    f.write('Visualize soma volume:'.ljust(lj) + '{}\n'.format(visualize_somavolumes))
    f.write('Root:'.ljust(lj) + '{}\n'.format(root))
    f.write('Graph folder:'.ljust(lj) + '{}\n'.format(infolder))
    f.write('Images folder:'.ljust(lj) + '{}\n'.format(infolder_tif))
    

#Write the outfile headers
with open(root + measurement_filename, 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow(['Condition'] + ['Age'] + ['Series'] + ['Name'] + ['Neurontype'] + 
                        ['Class'] + ['Length'] + ['MeanRadius'] + ['MaxRadius'] + ['Censored'])

with open(root + measurement_filename2, 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow(['Condition'] + ['Age'] + ['Series'] + ['Name'] + ['Neurontype'] + 
                        ['SomaOutgrowthLength'] + ['SomaOutgrowthCount'] + ['NeuriteOutgrowthLength'] + ['NeuriteOutgrowthCount'] + ['NeuriteOutgrowthPositionArray'] + ['MainBranchLength'] +
                        ['KinkCount'] + ['KinkPositionArray'] + ['SomaVolume'] + ['KinkDensityTotal'] + ['BeadCount'] + ['BeadDensity'] + ['PMVConnectionCount'] + ['BlobCount'] + ['L_R'] +
                        ['ThresholdAngle'] + ['Censored'])

        
for neurontype in neurontypes:
    n_infolder = root + neurontype + '/' + infolder
    n_infolder_tif = root + neurontype + '/' + infolder_tif
    n_outfolder_cleanedtrees = root + neurontype + '/' + outfolder_cleanedtrees
    n_outfolder_classifiedtrees = root + neurontype + '/' + outfolder_classifiedtrees
    n_outfolder_wavytrees = root + neurontype + '/' + outfolder_wavytrees
    if neurontype=='ALM' and visualize_somavolumes:
        n_outfolder_gifs = root + neurontype + '/' + outfolder_gifs
    
    
    #Make all the approriate folders
    if not os.path.exists(n_outfolder_cleanedtrees):
        os.makedirs(n_outfolder_cleanedtrees)
    if not os.path.exists(n_outfolder_classifiedtrees):
        os.makedirs(n_outfolder_classifiedtrees)
    if not os.path.exists(n_outfolder_wavytrees):
        os.makedirs(n_outfolder_wavytrees)
    try:
        if not os.path.exists(n_outfolder_gifs):
            os.makedirs(n_outfolder_gifs)
    except:
        n_outfolder_gifs=None
    
    files = os.listdir(n_infolder)
    files = [file for file in files if (file.endswith('.swc')) and file[0]!='#']
    for file in tqdm(files):
        file_tif = file[:-4] + '.tif'
        string = file[:-4]
        dat = string.split('_')
        strain = dat[0]
        series = dat[2]
        age = dat[1]
        name = dat[3]
        #gaussian_sigma = dat[4]
        #algorithm = dat[5]
        #decon = dat[6]
        
        raw_tree = utility.readSWC(n_infolder+file)
        
        if cleanup_tree == True:
            clean_tree = cleanup.cleanup(raw_tree,
                            neurontype=neurontype,
                            scale=scale,
                            visualize=False)
        else:
            clean_tree = raw_tree
            
        
        
        tree_lengths, tree_classes, tree_mean_radii, tree_max_radii, mainbranch, annotated_mainbranch, classified_tree = classify.classify(clean_tree,
                                                       neurontype=neurontype,
                                                       length_threshold=length_threshold,
                                                       plm_outgrowth_suppression=plm_outgrowth_suppression,
                                                       scale=scale)
        
        kinks_count, fully_annotated_mainbranch, kinks_swc, wavytree_swc = wavyness.wavyness(annotated_mainbranch,
                          angle_threshold=angle_threshold,
                          window_size_linear_regression=window_size_linear_regression,
                          window_size_maximum_supression=window_size_maximum_supression,
                          n_colors=10,
                          scale=scale,
                          fix_node=False,
                          plot_cdf=False)
        
        kink_positions, outgrowth_positions = kinkPositions.kinkPositions(fully_annotated_mainbranch, scale=scale)
        kink_positions_string = '-'.join([str(pos) for pos in kink_positions])
        outgrowth_positions_string = '-'.join([str(pos) for pos in outgrowth_positions])
        
        mainbranch_bead, bead_count = beads.countBeads(mainbranch, annotations=fully_annotated_mainbranch, scale=scale)
        
        
        density_total = kinks_count/tree_lengths[-1]
        
        

        #if False:
        if neurontype=='ALM':
            soma_volume = somaVolume.somavolume(classified_tree,
                                                filename_tif=n_infolder_tif+file_tif,
                                                scale=scale,
                                                visualize=visualize_somavolumes,
                                                outfolder_gif=n_outfolder_gifs,
                                                soma_volume_singleseed=soma_volume_singleseed)
        else:
            soma_volume = 0
        
        
        #Save results to .csv files
        with open(root + measurement_filename, 'a', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for i in range(len(tree_lengths)):
                spamwriter.writerow([strain] + [age] + [series] + [name] + [neurontype] + 
                                    [str(tree_classes[i])] + [str(tree_lengths[i])] + [str(tree_mean_radii[i])] + [str(tree_max_radii[i])]  + [str(0)])
                
        
        soma_total = 0
        neurite_total = 0
        for i in range(len(tree_lengths)):
            if tree_classes[i] == 'SomaOutgrowth':
                soma_total += tree_lengths[i]
            elif tree_classes[i] == 'NeuriteOutgrowth':
                neurite_total += tree_lengths[i]
        
        neurontype_plml_or_plmr = 'ALM' #ALM neurons can not be differentiated
        if neurontype=='PLM':
            if tree_classes.count('PMV')>0: #If a crossing PMV process was detected its a PLML neuron, else PLMR
                neurontype_plml_or_plmr = 'PLML'
            else:
                neurontype_plml_or_plmr = 'PLMR'

        
        with open(root + measurement_filename2, 'a', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
            mainbranch_length = utility.calculateDistancesTree(annotated_mainbranch, scale=scale)
            spamwriter.writerow([strain] + [age] + [series] + [name] + [neurontype] + 
                                [str(soma_total)] +  [str(tree_classes.count('SomaOutgrowth'))] +  [str(neurite_total)] + [str(tree_classes.count('NeuriteOutgrowth'))] + [outgrowth_positions_string]  +  [str(tree_lengths[-1])] +
                                [str(kinks_count)] + [kink_positions_string] + [str(soma_volume)] + [str(density_total)] + [str(bead_count)] + [str(bead_count/mainbranch_length)] + [str(tree_classes.count('PMV_Connection'))] +
                                [str(tree_classes.count('Blob'))] + [neurontype_plml_or_plmr] + [str(angle_threshold)] + [str(0)])
            #spamwriter.writerow([strain] + [age] + [name] + [series] + [neurontype] + ['SomaOutgrowth'] +  [str(soma_total)]  +  [str(tree_classes.count('SomaOutgrowth'))]  + [gaussian_sigma]  + [str(0)])        
        
        #save relevant .swc files
        utility.saveSWC(n_outfolder_classifiedtrees+file, classified_tree)
        utility.saveSWC(n_outfolder_wavytrees+file, wavytree_swc)
        utility.saveSWC(n_outfolder_wavytrees+file[:-4]+'_beads.swc', mainbranch_bead)
        #print('\n{}:\nBead count: {}\nKink count: {}\n\n'.format(file, int(bead_count), kinks_count))
                
MessageBeep()
end = time.time()
with open(root + parameters_filename, 'a') as f:
    lj = 30
    f.write('Total time [s]:'.ljust(lj) + '{:.2}\n'.format(end-start))