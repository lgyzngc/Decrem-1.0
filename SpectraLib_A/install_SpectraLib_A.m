% adds the asymmetric clustering and data generate/read directories to path
%
% $Authors: Marina Meila, Tatiana Maravina
% $Part of SpectraLib_A
% $Last revision: 06-June-2007
 
global SPECTRAL_HOME ASYMSPECTRAL_HOME
if isempty(ASYMSPECTRAL_HOME)
    ASYMSPECTRAL_HOME = '/costila/speclust/maravina/code-release';
end
addpath(genpath(ASYMSPECTRAL_HOME))

if isempty(SPECTRAL_HOME)
  disp('No directory listed for where Spectral Library stored.')
  disp('ASYMSPECTRAL needs Spectral library to function.  Please make sure it is installed before continuing.')
  SPECTRAL_HOME = '/costila/speclust/code-deepak-dec03';
end

addpath(genpath(SPECTRAL_HOME))
