function y=read_commomdata()
clear all; close all;

filename_contact = './contact10000.dat';
filename_coord = './coord10000.dat';
filename_allfaces = './allfaces10000.dat';
filename_list = './list10000.dat';

% filename_contact = './Lattice_protein_enumeration_code_and_structures/select-subset/contact103346.dat';
% filename_coord = './Lattice_protein_enumeration_code_and_structures/select-subset/coord103346.dat';
% filename_allfaces = './Lattice_protein_enumeration_code_and_structures/select-subset/allfaces103346.dat';


% filename_contact = './Lattice_protein_enumeration_code_and_structures/select-subset/contact10000.dat';
% filename_coord = './Lattice_protein_enumeration_code_and_structures/select-subset/coord10000.dat';
% filename_allfaces = './Lattice_protein_enumeration_code_and_structures/select-subset/allfaces10000.dat';
% 


A_coord=dlmread(filename_coord);

%plot coord
structure_ID = 8765;
%structure_ID = 873; 


x_array =  A_coord(structure_ID+1,1:3:end);
y_array =  A_coord(structure_ID+1,2:3:end);
z_array =  A_coord(structure_ID+1,3:3:end);
plot3(x_array,y_array,z_array,'-b');

A_coord(structure_ID+1,:),


for i=1:27
    text(x_array(i),y_array(i),z_array(i), num2str(i), 'Color', 'r','FontSize',20);
end

%stop; 



faces = zeros(6,9);
faces(1,:) = find(x_array==0); 
faces(2,:) = find(y_array==0);
faces(3,:) = find(z_array==0);
faces(4,:) = find(x_array==2); 
faces(5,:) = find(y_array==2);
faces(6,:) = find(z_array==2);
faces, 

A_allfaces=dlmread(filename_allfaces);
%construct AllFaces
AllFaces = A_allfaces;


curr_indexes = find(AllFaces(:,1)==structure_ID);
faces_fromfile = sort(AllFaces(curr_indexes(1:4:end),4:end)+1,2),
%faces_fromfile = sort(27-AllFaces(curr_indexes(1:4:end),4:end),2),





A_contact=dlmread(filename_contact);
curr_structure = structure_ID;
curr_indexes = find (A_contact(:,1)==curr_structure);
ContactMatrixLen(structure_ID)=length(curr_indexes);
ContactMatrixA{structure_ID} = A_contact(curr_indexes,2);
ContactMatrixB{structure_ID} = A_contact(curr_indexes,3);

%construct ContactMatrixLen ContactMatrixA and B
% unique_contact_structures = unique(A_contact(:,1));
% for unique_contact_structures_i=1:length(unique_contact_structures),
%     curr_structure = unique_contact_structures(unique_contact_structures_i);
%     curr_indexes = find (A_contact(:,1)==curr_structure);
%     ContactMatrixLen(unique_contact_structures_i)=length(curr_indexes);
%     ContactMatrixA{unique_contact_structures_i} = A_contact(curr_indexes,2);
%     ContactMatrixB{unique_contact_structures_i} = A_contact(curr_indexes,3);
% end

contacts_fromfile = [ContactMatrixA{structure_ID},ContactMatrixB{structure_ID}]+1,
%contacts_fromfile = 27-[ContactMatrixA{structure_ID},ContactMatrixB{structure_ID}],



stop;

A_list=dlmread(filename_list);





y=0;
end 
