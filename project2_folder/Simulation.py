import sys, os;
import numpy as np;
from mayavi import mlab;
import matplotlib.pyplot as mpl;



class Simulation:	
    def __init__(self,experiment_name,experiment_folder, config_file = None, N=8,makeCrystal='false', T_start=100, b=5.72, element="ar",T_bath=100, thermostatON='false', num_pro=4, max_time_steps=1000, num_time_steps=100, dt=0.01, writeVMD='false', writeMeasurements='true', writeToFilename = None, Berendsen='false', Andersen='false'):
        self.experiment_name = experiment_name; #typical make_run_thermalize
        self.experiment_folder = experiment_folder;
        self.measurement_file_name = experiment_name;
        self.current_part = 0;
        self.measure_spatial_pressure_at = -1;
        self.num_thermostatON = 0;
        self.num_solve = num_time_steps;
        self.measurement_frequency = 1;
        self.Fx = 0;
        self.total_num_time = 0;
        
 
        self.state_file = experiment_folder + experiment_name+'_state%d.xyz'%self.current_part;
        self.program_folder = '/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/Project2-build-Desktop-Release/'
        self.executable = self.program_folder+'Project2'
        self.max_time_steps = max_time_steps;
        self.real_current_time_step = 0;
        if(config_file==None):
            self.config_file = experiment_folder + experiment_name + '_config%d.cfg'%self.current_part;
            self.read_attributes(N,makeCrystal, T_start, b, element,T_bath, thermostatON, num_pro, max_time_steps, num_time_steps, dt, writeVMD, writeMeasurements, writeToFilename, Berendsen, Andersen);
        else:
            self.config_file = config_file;
            self.read_config_file();
                
    def read_attributes(self,N,makeCrystal, T_start, b, element,T_bath, thermostatON, num_pro, max_time_steps, num_time_steps, dt, writeVMD, writeMeasurements, writeToFilename, Berendsen, Andersen):
        L0 = 3.405;
        self.N = N;
        self.makeCrystal=makeCrystal;
        self.T_start = T_start;
        self.b=b;
        self.element = element;
        self.T_bath = T_bath;
        self.thermostatON = thermostatON;
        self.num_pro = num_pro;
        self.max_time_steps = max_time_steps;
        self.num_time_steps = num_time_steps;
        self.dt = dt;
        self.writeVMD = writeVMD;
        self.writeMeasurements = writeMeasurements;
        self.writeToFilename = writeToFilename;
        self.Berendsen = Berendsen;
        self.Andersen = Andersen;
        r_cut = 3*L0;
        self.CellN = int(N*b/r_cut);
    	self.r_cut = N*b/self.CellN;
        
    def read_config_file(self):
        infile = open(self.config_file, 'r');
        lines = infile.readlines();
        for i in range(len(lines)):
            line = lines[i];
            if(line[0]!='/'):
                line=line.replace("true", "'true'")
                line=line.replace("false", "'false'")
                exec "self."+line
            
    def set_state_file(self, newfile):
        new_state_file = open(newfile, 'r');
        filestring = new_state_file.read();
        this_state_file = open(self.state_file, 'w');
        this_state_file.write(filestring)
        #self.firstStateFile = newfile;
    def set_N(self, N):
        L0 = 3.405;
        self.N = N;
        r_cut = 3*L0;
        self.CellN = int(N*self.b/r_cut);
    	self.r_cut = N*self.b/self.CellN;
        
    def generate_config_file(self):
        next_state_num = self.current_part +1;
        next_ending = '_state%d.xyz";\n'%next_state_num
        print "current part", self.current_part;
        self.config_file = self.experiment_folder + self.experiment_name + '_config%d.cfg'%self.current_part;
        outfile = open(self.config_file, 'w')
        outfile.write('//stupid comment\n');
        outfile.write('makeCrystal = %s;\n'%self.makeCrystal)
        #outfile.write('firstStateFile = "'+self.firstStateFile+'";\n')
        outfile.write('state_file = "' + self.state_file + '";\n');
        outfile.write('next_state_file = "'+ self.experiment_folder + self.experiment_name + next_ending)
        outfile.write('N = %d;\n'%self.N)
        outfile.write('measure_spatial_pressure_at = %d;\n'%self.measure_spatial_pressure_at);
        outfile.write('num_thermostatON = %d;\n'%self.num_thermostatON)
        outfile.write('CellN = %d;\n'%self.CellN)
        outfile.write('r_cut = %f;\n'%self.r_cut)
        outfile.write('T_start = %f;\n'%self.T_start)
        outfile.write('T_bath = %f;\n'%self.T_bath)
        outfile.write('b = %f;\n'%self.b)
        outfile.write('Fx = %f;\n'%self.Fx);
        outfile.write('thermostatON = %s;\n'%self.thermostatON)
        outfile.write('element = "%s";\n'%self.element)
        outfile.write('num_pro = %d;\n'%self.num_pro);
        outfile.write('max_time_steps = %d;\n'%self.max_time_steps);
        outfile.write('num_time_steps = %d;\n'%self.num_time_steps)
        outfile.write('dt = %f;\n'%self.dt)
        outfile.write('writeVMD = %s;\n'%self.writeVMD);
        outfile.write('writeMeasurements = %s;\n'%self.writeMeasurements);
        outfile.write('writeToFilename = "%s";\n'%self.writeToFilename);
        outfile.write('Berendsen = %s;\n'%self.Berendsen);
        outfile.write('Andersen = %s;\n'%self.Andersen);
        outfile.write('measurement_file_name = "%s";\n'%self.measurement_file_name);
        outfile.write('measurement_frequency = %d;\n'%self.measurement_frequency);
        outfile.write('real_current_time_step = %d;\n'%self.real_current_time_step);
        outfile.close();
    	
    def make_crystal(self,N, T_start,b,element):
        L0 = 3.405;
    	r_cut = 3*L0;
        self.N=N;
        self.T_start=T_start;
        self.b=b;
        self.element = element;
    	self.CellN = int(N*b/r_cut);
    	self.r_cut = N*b/self.CellN;
        self.makeCrystal = 'true';
        self.generate_config_file();
        self.read_config_file();
        #self.run();
        	
    def calculate_permeability(self,U,mu,phi,n):
        L0 = 3.405;
        self.a= (U*mu/(n*self.Fx));
        self.a = self.a*(1e-10*L0)**2
        return self.a
		
    def thermalize(self, T, N, writeVMD = None):
        #run the program with thermostat for N timesteps
        if(writeVMD!=None):
            self.writeVMD = writeVMD;
        self.T_bath = T;
        self.num_thermostatON = N;
        self.thermostatON= 'true';
        self.Berendsen = 'true';
        #self.generate_config_file();
        #self.makeCrystal = 'false';
        #self.run();
        #self.thermostatON='false'
    def read_state_file(self):
        state_file = open(self.start_state_file, 'r')
        lines = state_file.readlines();
        positions = [];
        velocities = [];
        for i in range(2,len(lines)):
            line = lines[i];
            words = line.split();
            x = eval(words[1]);
            y = eval(words[2]);
            z = eval(words[3]);
            vx = eval(words[4]);
            vy = eval(words[5]);
            vz = eval(words[6]);
            positions.append(array([x,y,z]));
            velocities.append(array([vx,vy,vz]));
        return positions,velocites;

        
    def makePores(self, num, r_min, r_max):
        #make random spheres of random size
        position_list = [];
        radius_list = np.zeros(num);
        L0 = 3.405;
        r_min = r_min/(L0);
        r_max = r_max/(L0);
        L = self.N*self.b/L0;
        print L
        #find the pores
        for i in range(num):
            position_list.append(np.random.uniform(0,L,3));
            radius_list[i] = np.random.uniform(r_min, r_max,1);
        #run through the xyz-file and change canMove
        #self.state_file = self.experiment_folder + self.experiment_name + '_config%d.cfg'%self.current_part;
        self.state_file = self.experiment_folder + self.experiment_name+'_state%d.xyz'%self.current_part;
        state_file = open(self.state_file, 'r')
        lines = state_file.readlines();
        state_file.close()
        for i in range(2,len(lines)-1):
            line = lines[i];
            words = line.split();
            x = eval(words[1]);
            y = eval(words[2]);
            z = eval(words[3]);
            for pos, r in zip(position_list, radius_list):
                R_sq = (pos[0]-x)**2 + (pos[1]-y)**2  + (pos[2]-z)**2;
                if(R_sq< r**2):
                    words[-1] = "0\n";
                    lines[i] = " ".join(words);
                    break
                
        new_state_file = open(self.state_file, 'w');
        filestring = "".join(lines);
        new_state_file.write(filestring);
        print "the pores are in ", self.state_file
    def makeCyllinder(self, rel_pos_x, rel_pos_y,radius):
        L0 = 3.405;
        R = radius/L0;
        state_file = open(self.state_file, 'r');
        lines = state_file.readlines();
        state_file.close();
        L = self.N*self.b/L0;
        x_center = rel_pos_x*L;
        y_center = rel_pos_y*L;
        for i in range(2, len(lines)-1):
            line = lines[i];
            words = line.split();
            x = eval(words[1]);
            y = eval(words[2]);
            z = eval(words[3]);
            radius_sq = (x-x_center)**2 + (y-y_center)**2;
            if(radius_sq>R**2):
                words[-1] = "0\n";
                lines[i] = " ".join(words);
        new_state_file = open(self.state_file, 'w');
        filestring = "".join(lines);
        new_state_file.write(filestring);
        print "made cyllinder in ", self.state_file
        
    def set_external_force(self, Fx):
        self.Fx = Fx;
    def run(self):
        self.measurement_file_name = self.experiment_name + "%d"%self.current_part
        self.total_num_time +=self.num_time_steps;
        self.generate_config_file();
        self.current_part +=1;
        print "run using " + self.config_file+ " with thermo:"+ self.thermostatON
        print "reading from ", self.state_file;
        #run the program with the command_line_file.
        os.system(self.executable + ' '+ self.experiment_name + " " + self.experiment_folder + " " + self.config_file)
        self.real_current_time_step +=self.num_time_steps;
        #self.state_file = self.experiment_folder + self.experiment_name + "_state%d.xyz"%self.current_part
        self.state_file = self.experiment_folder + self.experiment_name + "_state%d.xyz"%self.current_part

    def solve(self, N, writeVMD=None):
        #run the program with N timesteps
        #the xyz files get filename_ending solve to separate from files made during thermalize
        if(writeVMD!=None):
            self.writeVMD=writeVMD;
        self.num_time_steps = N;
        #self.generate_config_file();
        self.run();
    def glueVMDfiles(self):
        #makes two files, one glued thermalize..., and one glued from running solve
        glue_folder = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/Simulations/"
        glue_file_name = "glueVMDnew.py"
        os.system("python " + glue_folder+glue_file_name + " " + self.experiment_name + " "+ self.experiment_folder + " %d rm"%(self.real_current_time_step));
    def glueMeasurements(self, num_parts):
        glue_folder = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/Simulations/";
        glue_file_name = "glueMeasurements.py"
        os.system("python " + glue_folder+glue_file_name + " "+ self.experiment_name + " " + self.experiment_folder + " %d rm"%num_parts);
    def reduce_fluid_density_to(self, fraction):
        #run through the xyz-file and remove half the particles with canMove=1
        state_file = open(self.state_file, 'r')
        lines = state_file.readlines();
        new_lines = [];
        new_lines.append(lines[0]);
        new_lines.append(lines[1]);
        j=0;
        state_file.close()
        for i in range(2,len(lines)-1):
            append = True;
            line = lines[i];
            words = line.split();
            if eval(words[-1])==1:
                draw = np.random.random();
                if(draw>fraction):
                    append=False;
            if(append):
                j+=1;
                new_lines.append(lines[i]);
        new_lines[0] = "%d\n"%j;
        new_state_file = open(self.state_file, 'w');
        filestring = "".join(new_lines);
        print "j: ",j;
        new_state_file.write(filestring);
        
    def compile_program(self):
        #compiles. Should only happen once
        terminal_compile = 'qmake-qt4 /home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/Project2/Project2.pro -r -spec linux-g++-64'
        terminal_compile2 = 'make -w /home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/Project2-build-Desktop-Release'
        os.system(terminal_compile);
        os.system(terminal_compile2);
    def plot_measurements(self):
        plot_program = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/results/plot_energy.py"
        os.system("python " + plot_program + " glued_" + self.experiment_name + " " + self.experiment_folder + "/results/  all 1");
    def plot_spatial_pressure(self,layer):
        spatial_file_name = self.experiment_folder + "/results/" + self.experiment_name + "_spatial_pressure.txt";
        spatial_file = open(spatial_file_name, 'r');
        file_lines = spatial_file.readlines();
        surf = np.zeros([self.CellN, self.CellN]);
        index1 = layer*self.CellN*self.CellN;
        for i in range(index1, index1 + self.CellN):
            words = file_lines[i].split();
            for j in range(self.CellN):
                surf[i][j] = eval(words[j]);
        X,Y = np.meshgrid(np.linspace(0,self.CellN-1, self.CellN), np.linspace(0,self.CellN-1, self.CellN))
        s = mlab.mesh(X,Y, surf);
        mlab.colorbar(s, orientation='vertical');
        axem = mlab.axes();
        outlinem = mlab.outline();
        mlab.xlabel('x');
        mlab.ylabel('y');
        mlab.zlabel('z');
        mlab.show();
    def plot_flow_profile(self, start_time_step, stop_time_step, rel_centerX, rel_centerY, max_radius, num_bins=100):
        result_folder = self.experiment_folder + "/results/";
        L0 = 3.405 #aangstrom
        t0 = 2.1569e3#fs
        bin_size = max_radius/(num_bins*L0);
        print "bin size", bin_size
        distribution = np.zeros(num_bins);
        flow_profile = np.zeros(num_bins);
        density = np.zeros(num_bins);
        r_plot = np.linspace(0, max_radius, num_bins);
        glued_file_str = result_folder + "glued_" + self.experiment_name + ".xyz";
        print "reading glued file"
        #glued_file = open(glued_file_str, 'r');
        #all_lines = glued_file.readlines();
        #glued_file.close();
        L = self.N*self.b/L0;
        centerX = rel_centerX*L;
        centerY = rel_centerY*L;
        start_line = 0;
        stop_line = 0;
        #dl = (len(all_lines))/float(self.total_num_time);
        #print dl
        for j in range(start_time_step, stop_time_step):
            current_state_file_name = self.experiment_name + "%d.xyz"%j;
            current_file = open(result_folder + current_state_file_name, 'r');
            lines = current_file.readlines();
            current_file.close();
            #start = int (j*dl);
            #stop = int ((j+1)*dl);
            #lines = all_lines[start:stop];
            num = 0;
            for i in range(2,len(lines)-1):
                words = lines[i].split();
                canMove = eval(words[-1]);
                if(canMove):
                    x = eval(words[1]);
                    y = eval(words[2]);
                    z = eval(words[3]);
                    vz = eval(words[6]);
                    r = np.sqrt((x-centerX)**2 + (y-centerY)**2)
                    if(r<max_radius/L0):
                        #  shell_volume = 4*np.pi*r**2*bin_size;
                        distribution[int(r/bin_size)]+=1;
                        #density[int(r/bin_size)] +=1/shell_volume;
                        flow_profile[int(r/bin_size)]+=vz;
                        #print "doing stuff", j
                    num+=1;
            print num;
        volumes = (((r_plot+bin_size))**2 - (r_plot)**2)*L/L0**2*np.pi;
        total_num = sum(distribution);
        total_volume = L*np.pi*(max_radius/L0)**2;
        viscosities = np.zeros(num_bins)
        a = 20/L0 #aa
        F = self.Fx#*3.0303*10**(-1)
        self.n = total_num/total_volume;
        print n,F
        for i in range(num_bins):
            flow_profile[i] = (flow_profile[i]/float(distribution[i]));
            density[i] = distribution[i]/(volumes[i]);
            viscosities[i] = (a**2-(r_plot[i]/L0)**2)/(4*flow_profile[i])*n*F
        viscosity = sum(viscosities)/len(viscosities);
        print "viscosity = ", viscosity, "md"#"eVns/AA^3"

        #radii = np.linspace()
        #density = distribution/volumes;
        mpl.plot(r_plot, flow_profile)
        mpl.title('Flow profile')
        mpl.xlabel('radius (AAngstrom)')
        mpl.ylabel('velocity (AA/ns)');
        mpl.figure()
        mpl.plot(r_plot, density);
        mpl.xlabel('radius (AAngstrom)')
        mpl.ylabel('number density (AA^-3)')
        mpl.title('Radial number density');
        mpl.show();
            
            


                
            

if __name__=='__main__':
    simulator = Simulation(6);
    simulator.compile_program();
    simulator.run();

