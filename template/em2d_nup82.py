"""@namespace IMP.pmi.restraints.em2d
Restraints for handling electron microscopy images.
"""

from __future__ import print_function
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.em2d
import IMP.pmi.tools
import ihm.location
import ihm.dataset


#############################################################################################
# PCAFitRestraint - compares how well the principal components of the segmented class average
# fit to the principal components of the particles, developed by Dina Schneidman-Duhovny.
#############################################################################################
class ElectronMicroscopy2D():

    def __init__(
        self,
        representation,
        images,
        pixel_size,
        image_resolution,
        projection_number,
        micrographs_number=None,
        resolution=None,
        n_components = 1):

        self.weight=1.0
        self.m = representation.prot.get_model()
        self.rs = IMP.RestraintSet(self.m, 'em2d')
        self.label = "None"

        self.datasets = []
        for image in images:
            if representation:
                d = representation.get_file_dataset(image)
                if d:
                    self.datasets.append(d)
                    continue
            l = ihm.location.InputFileLocation(image,
                                 details="Electron microscopy class average")
            d = ihm.dataset.EM2DClassDataset(l)
            self.datasets.append(d)

        if representation:
            for p, state in representation._protocol_output:
                for i in range(len(self.datasets)):
                    p.add_em2d_restraint(state, self, i, resolution, pixel_size,
                                         image_resolution, projection_number,
                                         micrographs_number)

        ps = []

        ################################
        # Dyn2.1
        ################################
        particles = IMP.pmi.tools.select(representation,
                                        resolution=resolution,
                                        name="Dyn2.1")
        for p in particles:
            ps.append(p)

        # Dyn2.2
        particles = IMP.pmi.tools.select(representation,
                                        resolution=resolution,
                                        name="Dyn2.2")
        for p in particles:
            ps.append(p)


        ################################
        # Nup82.1
        ################################
        particles = IMP.pmi.tools.select(representation,
                                        resolution=resolution,
                                        name="Nup82.1")
        for p in particles:
            ps.append(p)

        # Nup82.2
        particles = IMP.pmi.tools.select(representation,
                                        resolution=resolution,
                                        name="Nup82.2")
        for p in particles:
            ps.append(p)


        ################################
        # Nup159.1
        ################################
        particles = IMP.pmi.tools.select(representation,
                                        resolution=resolution,
                                        name="Nup159.1")
        for p in particles:
            self_names=(p.get_name()).replace("-","_").split("_")

            if int(self_names[1]) > 1116 :
                ps.append(p)

        # Nup159.2
        particles = IMP.pmi.tools.select(representation,
                                        resolution=resolution,
                                        name="Nup159.2")
        for p in particles:
            self_names=(p.get_name()).replace("-","_").split("_")

            if int(self_names[1]) > 1116 :
                ps.append(p)


        ################################
        # Nsp1.1
        ################################
        particles = IMP.pmi.tools.select(representation,
                                        resolution=resolution,
                                        name="Nsp1.1")
        for p in particles:
            #parent_names=(p.get_parent().get_name()).replace("-","_").split("_")
            self_names=(p.get_name()).replace("-","_").split("_")

            if int(self_names[1]) > 636 :
                ps.append(p)

        # Nsp1.2
        particles = IMP.pmi.tools.select(representation,
                                        resolution=resolution,
                                        name="Nsp1.2")
        for p in particles:
            self_names=(p.get_name()).replace("-","_").split("_")

            if int(self_names[1]) > 636 :
                ps.append(p)

        print ("len(ps) = ", len(ps))   #len(particles_all)=3174    #len(particles_selected)=2076
        particles = ps
        ################################
        
        map = IMP.em.SampledDensityMap(ps, image_resolution, pixel_size)
        IMP.em.write_map(map, "em2d_particle_selected.mrc")

        """
        # IMP.atom.get_by_type
        particles = IMP.pmi.tools.select(
            representation,
            resolution=resolution)
        #print particles
        """
        # read PGM FORMAT images
        # format conversion recommendataion - first run "e2proc2d.py $FILE ${NEW_FILE}.pgm"
        # then, run "convert ${NEW_FILE}.pgm -compress none ${NEW_FILE2}.pgm"
        if (n_components >= 2) :    # Number of the largest components to be considered for the EM image
            em2d = IMP.em2d.PCAFitRestraint(
                particles, images, pixel_size, image_resolution, projection_number, True, n_components)
        else :
            em2d = IMP.em2d.PCAFitRestraint(
                particles, images, pixel_size, image_resolution, projection_number, True)
        self.rs.add_restraint(em2d)

    def set_label(self, label):
        self.label = label

    def add_to_model(self):
        IMP.pmi.tools.add_restraint_to_model(self.m, self.rs)

    def get_restraint(self):
        return self.rs

    def set_weight(self,weight):
        self.weight=weight
        self.rs.set_weight(self.weight)

    def get_output(self):
        self.m.update()
        output = {}
        score = self.weight*self.rs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["ElectronMicroscopy2D_" + self.label] = str(score)
        return output


#############################################################################################
# EM2DRestraint_FFT - implements FFT based image alignment, developed by Javier Velazquez-Muriel
#############################################################################################
class ElectronMicroscopy2D_FFT():

    def __init__(
        self,
        representation,
        images,
        pixel_size,
        image_resolution,
        projection_number,
        resolution=None):

        self.weight=1.0
        self.m = representation.prot.get_model()
        self.rs = IMP.RestraintSet(self.m, 'em2d_FFT')
        self.label = "None"

        # read SPIDER FORMAT images (format conversion should be done through EM2EM)
        srw = IMP.em2d.SpiderImageReaderWriter()
        imgs = IMP.em2d.read_images(images, srw)
        rows = imgs[0].get_header().get_number_of_rows()
        cols = imgs[0].get_header().get_number_of_columns()

        # pixel_size: sampling rate of the available EM images
        # image_resolution: resolution at which you want to generate the projections of the model
        #                   In principle you want "perfect" projections, so use the highest resolution
        # projection_number: Number of projections to use for the initial registration
        #                    (coarse registration) to estimate the registration parameters
        params = IMP.em2d.Em2DRestraintParameters(pixel_size, image_resolution, projection_number)

        # This method (recommended) uses preprocessing of the images and projections to speed-up the registration
        params.coarse_registration_method = IMP.em2d.ALIGN2D_PREPROCESSING
        params.optimization_steps = 50
        params.simplex_initial_length = 0.1
        params.simplex_minimum_size = 0.02

        # use true if you want to save the projections from the model that best match the EM images
        params.save_match_images = False

        ######################
        # set up the em2D restraint
        ######################
        score_function = IMP.em2d.EM2DScore()
        em2d_restraint = IMP.em2d.Em2DRestraint(self.m)
        em2d_restraint.setup(score_function, params)
        em2d_restraint.set_images(imgs)
        em2d_restraint.set_fast_mode(5)
        em2d_restraint.set_name("em2d_restraint")

        # IMP.atom.get_by_type
        particles = IMP.pmi.tools.select(
            representation,
            resolution=resolution)
        print ("len(particles) = ", len(particles))
        """
        ps = []
        for p in particles:
            ps.append(p)
            print p
            print IMP.core.XYZR(p).get_radius()
            print IMP.atom.Mass(p).get_mass()
        """
        container = IMP.container.ListSingletonContainer(self.m, particles)
        em2d_restraint.set_particles(container)

        self.rs.add_restraint(em2d_restraint)

    def set_label(self, label):
        self.label = label

    def add_to_model(self):
        IMP.pmi.tools.add_restraint_to_model(self.m, self.rs)

    def get_restraint(self):
        return self.rs

    def set_weight(self,weight):
        self.weight=weight
        self.rs.set_weight(self.weight)

    def get_output(self):
        self.m.update()
        output = {}
        score = self.weight*self.rs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["ElectronMicroscopy2D_FFT_" + self.label] = str(score)
        return output
