#!/usr/bin/env python

from __future__ import print_function
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.pmi1.restraints.crosslinking
import IMP.pmi1.restraints.stereochemistry
import IMP.pmi1.restraints.em
import IMP.pmi1.restraints.em2d
import IMP.pmi1.restraints.basic
import IMP.pmi1.restraints.proteomics
import IMP.pmi1.representation
import IMP.pmi1.tools
import IMP.pmi1.samplers
import IMP.pmi1.output
import IMP.pmi1.macros
import os
import sys

from math import sqrt
import IMP.em2d as em2d
import IMP.rmf
import RMF
import logging
import linecache
import time

####sys.settrace

def score_model(rmffile, image_file, pixel_size, n_projections, resolution, image_resolution,frame_number):
    IMP.base.set_check_level(IMP.base.USAGE_AND_INTERNAL)
    
    images = []
    with open(image_file, "r") as ins:
        for line in ins:
            images.append(str(line).strip('\n'))

    srw = em2d.SpiderImageReaderWriter()
    imgs =  em2d.read_images(images, srw)
    rows = imgs[0].get_header().get_number_of_rows()
    cols = imgs[0].get_header().get_number_of_columns()

    m = IMP.kernel.Model()
    inf = RMF.open_rmf_file_read_only(rmffile)
    h = IMP.rmf.create_hierarchies(inf, m)[0]
    particle2s = IMP.core.get_leaves(h)
    IMP.rmf.load_frame(inf, frame_number)
    children = h.get_children()

    ps = []
    for p in particle2s:
        if p.get_parent().get_name() == "Beads" and p.get_name()[0:6] == "Dyn2.1": # and p.get_parent().get_name()[0:6] == "Dyn2.1":
            ps.append(p)
        elif p.get_parent().get_name()[0:6] == "Dyn2.1" and p.get_name()[0:1] != "D":
            ps.append(p)
        elif p.get_parent().get_name() == "Beads" and p.get_name()[0:6] == "Dyn2.2":
            ps.append(p)
        elif p.get_parent().get_name()[0:6] == "Dyn2.2" and p.get_name()[0:1] != "D":
            ps.append(p)
        elif p.get_parent().get_name() == "Beads" and p.get_name()[0:7] == "Nup82.1":
            ps.append(p)
        elif p.get_parent().get_name()[0:7] == "Nup82.1" and p.get_name()[0:1] != "N":
            ps.append(p)
        elif p.get_parent().get_name() == "Beads" and p.get_name()[0:7] == "Nup82.2":
            ps.append(p)
        elif p.get_parent().get_name()[0:7] == "Nup82.2" and p.get_name()[0:1] != "N":
            ps.append(p)
        elif p.get_parent().get_name() == "Beads" and p.get_name()[0:10] == "Nup159.1_1" and int(p.get_name()[9:13]) > 1116:
            ps.append(p)
        elif p.get_parent().get_name()[0:10] == "Nup159.1_1"  and p.get_name()[0:1] != "N":
            ps.append(p)
        elif p.get_parent().get_name() == "Beads" and p.get_name()[0:10] == "Nup159.2_1" and int(p.get_name()[9:13]) > 1116:
            ps.append(p)
        elif p.get_parent().get_name()[0:10] == "Nup159.2_1"  and p.get_name()[0:1] != "N":
            ps.append(p)
        elif p.get_parent().get_name() == "Beads" and p.get_name()[0:8] == "Nsp1.1_7":
            ps.append(p)
        elif p.get_parent().get_name()[0:6] == "Nsp1.1" and p.get_name()[0:1] != "N":
            ps.append(p)
        elif p.get_parent().get_name() == "Beads" and p.get_name()[0:8] == "Nsp1.2_7":
            ps.append(p)
        elif p.get_parent().get_name()[0:6] == "Nsp1.2" and p.get_name()[0:1] != "N":
            ps.append(p)
        elif p.get_parent().get_name() == "Beads" and p.get_name()[0:3] == "662":
            ps.append(p)

    #IMP.rmf.add_hierarchies(ouf, ps)
    #IMP.rmf.save_frame(ouf,str(0))
    #map = IMP.em.SampledDensityMap(ps, 20, 3.225)
    #IMP.em.write_map(map, "Example.mrc")
    

    proj_params = em2d.get_evenly_distributed_registration_results(n_projections)
    opts = em2d.ProjectingOptions(pixel_size, resolution)
    projections = em2d.get_projections(ps,proj_params,rows,cols,opts)
    if len(projections) != n_projections:
        print("Problem generating projections")


    finder = em2d.ProjectionFinder()
    score_function = em2d.EM2DScore()
    subjects = em2d.read_images(images, srw)

    OPT_STEPS = int(100)
    INITIAL_LENGTH = 0.1
    MINIMUM_SIZE = 0.01

    #print("We are refining the projections!")
    params = em2d.Em2DRestraintParameters(pixel_size, resolution, n_projections)
    #print(n_projections)
    params.coarse_registration_method = em2d.ALIGN2D_PREPROCESSING
    params.optimization_steps = OPT_STEPS
    params.simplex_initial_length = INITIAL_LENGTH
    params.simplex_minimum_size = MINIMUM_SIZE
    params.save_match_images = True

    init_time = time.time()
    finder.setup(score_function, params)
    finder.set_model_particles(ps)
    finder.set_subjects(subjects)
    finder.set_projections(projections)
    finder.set_fast_mode(2)
    finder.get_complete_registration()

    #print("We are done registering! We now need to restore the scores!")
    all_registration_results = []
    registration_results = finder.get_registration_results()
    regnumber=0
    for reg in registration_results:
        #print(reg, regnumber)
        all_registration_results.append(reg)
        regnumber+=1

    for i in range(0, len(registration_results)):
        imgx = em2d.Image()
        imgx.set_size(rows, cols)
        em2d.get_projection(imgx, ps,registration_results[i], opts)
        ccc = em2d.get_cross_correlation_coefficient(subjects[i].get_data(),imgx.get_data())
        print(rmffile, i, "ccc=", ccc)

    #print("We are done here! I hope you had fun waiting!")
    em2d.write_registration_results("Registration-Parameters", all_registration_results)
    timelasted = time.time() - init_time
    #print("score_model: time complete registration %s" % timelasted)
    #print("coarse registration time %s", finder.get_coarse_registration_time())
    #print("fine registration time %s", finder.get_fine_registration_time())
    return all_registration_results
    
    
def parse_args():
    parser = IMP.OptionParser(usage="""%prog RMF Image-File PixelSize NumberOfProjections Resolution Image_Resolution
RMF: RMF files for loading the frames 
ImageFile: Selection file containing the images used for scoring
PixelSize:  Pixel size of the images in Angstrom/pixel
NumberOfProjection:   Number of projections of the model used for start registering
Resolution:     Resolution to generate the projections
Image_Resolution:    Images per batch used when scoring a lot of images (to avoid memory problems)
Frame Number : Frame nuimber to do . This is especially useful on the cluster for hihg-thourhgput score. 
""",
                              description="Filtering the models based on 2D Class Averages",
                              imp_module=em2d)
    opts, args = parser.parse_args()
    if len(args) != 7:
        parser.error("Wrong number of arguments! Check arguments lists.")
    return args

if __name__ == "__main__":
    #print("Program is starting. Do not worry, It will be over soon!\n")
    IMP.base.set_log_level(IMP.base.VERBOSE)
    args = parse_args()
    rmffile = args[0]
    image_file = args[1]
    pixel_size = float(args[2])
    n_projections = int(args[3])
    resolution = float(args[4])
    image_resolution = float(args[5])
    frame_number = int(args[6])
    
    all_regs = score_model(rmffile, image_file, pixel_size, n_projections, resolution, image_resolution, frame_number)
    
    #for i, reg in enumerate(all_regs):
    #print("Score for image %s: %f" % (i, reg.get_score()))
    #print("%s Global score %s" % (rmffile, em2d.get_global_score(all_regs)))
