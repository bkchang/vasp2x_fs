%vasp2x_fs

% This script is used for transforming VASP output into .bxsf of the 
% XCrysDen Fermi surface format.
% 
% One should first use VASP with ISYM=0, IBRION = 0, and a dense
% Gamma-centered Monkhorst-Pack grid to conduct a non-scf calculation, 
% then put the following files into the same directory and run the script.
%
% Files needed:
% EIGENVAL and POSCAR from your VASP output
% import_eigenval.m, import_lattice.m, reciprocal_lattice.m from the
% VASPLAB package, which could be found at
% https://www.mathworks.com/matlabcentral/fileexchange/36836-vasplab
%
% One should set the filename, kpoints, and the Fermi energy at the 
% parameter setting section below. Once the job is done, a file with the
% name "${filename}.bxsf" will be generated under the same directory. 
% One can then run "xcrysden --bxsf ${filename}.bxsf" to see the Fermi
% surface plot in XCrysDen.
%
% Author: Benjamin K. Chang,
% Institute of Atomic and Molecular Sciences, Academia Sinica, Taiwan
% email address: bkchang8@gmail.com
% Website: https://github.com/bkchang
% July 2017; Last revision: 27-July-2017

%------------- BEGIN CODE --------------

%% Setting Parameters
filename = 'name';
num_k_x = 24;
num_k_y = 24;
num_k_z = 24;
Ef = 6.000;       % Fermi energy in eV

%% Read from EIGENVAL
[ eigenvalues, kpoints, ~ ] = import_eigenval('EIGENVAL');
[ num_k, num_band ] = size(eigenvalues);

%% Read from POSCAR
POSCAR = import_poscar('POSCAR');
lattice = [0 0 0; reciprocal_lattice(POSCAR.lattice)];

%% Check the input kpoints
if num_k ~= num_k_x*num_k_y*num_k_z
    disp('The given kpoints are not consistent with the EIGENVAL file. Exiting script.');
    return
end

%% Combine the kpoints and the eigenvalues
k_eigen = [ kpoints(:,1:3) eigenvalues ];


%% Normalize the kpoints into the first quadrant in the k space (the
% reciprocal lattice cell)
for i=1:num_k
    if k_eigen(i,1) < 0
        k_eigen(i,1) = k_eigen(i,1) + 1;
    end
    if k_eigen(i,2) < 0
        k_eigen(i,2) = k_eigen(i,2) + 1;
    end
    if k_eigen(i,3) < 0
        k_eigen(i,3) = k_eigen(i,3) + 1;
    end
end

%% Sort the data by the k_z, k_y, k_x values, subsequently
k_eigen = sortrows(k_eigen, [1 2 3]);   
eigen = k_eigen(:, 3+1:3+num_band);

%% Save the blockwise arranged data into the variable "data"
data = [];
start_datarow = [];
element_per_layer = num_k_y * num_k_z;
for i_band = 1:num_band
    for i_layer = 1:num_k_x
        start_datarow = element_per_layer * (i_layer-1);
        data(:, :, i_layer, i_band) = reshape( eigen( start_datarow+1:start_datarow+element_per_layer, i_band), [ num_k_z num_k_y ])';
    end
end

%% Writing file
fid = fopen(strcat(filename, '.bxsf'), 'at');
fprintf(fid,'BEGIN_INFO\r\n');
fprintf(fid, '  Fermi Energy: %.4f \r\n', Ef);
fprintf(fid,'END_INFO\r\n');
fprintf(fid,'BEGIN_BLOCK_BANDGRID_3D\r\n');
fprintf(fid, strcat(filename, '\r\n'));
fprintf(fid,'BEGIN_BANDGRID_3D\r\n');
fprintf(fid, '    %d\r\n', num_band);
fprintf(fid, '    %d %d %d\r\n', num_k_x, num_k_y, num_k_z);
dlmwrite(strcat(filename,'.bxsf'), lattice,'delimiter',' ','precision',4,'-append');
fprintf(fid,'\r\n');
for ind1 = 1:num_band
    fprintf(fid,'   BAND:   %d\r\n', ind1);
    for ind2 = 1:num_k_x
        dlmwrite(strcat(filename,'.bxsf'), data(:, :, ind2, ind1),'delimiter',' ','precision',4,'-append')
        fprintf(fid,'\r\n');
    end
end
fprintf(fid,'END_BANDGRID_3D\r\n');
fprintf(fid,'END_BLOCK_BANDGRID_3D');
fclose(fid);
disp('Job done.')
%------------- END OF CODE --------------