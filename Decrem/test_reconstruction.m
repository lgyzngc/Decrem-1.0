% test to generate decoupled metabolic model for bacillus iYO844 model
model = load('..\bacillus\iYO844.mat','-mat');
cofactor_path = '..\bacillus\cofactor.txt';
secrated_path = '..\bacillus\secrated_bacillus.txt';
nutrient_path = '..\bacillus\nutrient_bacillus.txt'
general_IO_path = '..\bacillus\general_IO_bacillus.txt'

reconstructed_model = decoupledModelConstruct(model,cofactor_path,secrated_path,nutrient_path,general_IO_path)