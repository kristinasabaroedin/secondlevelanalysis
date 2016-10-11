clear all; close all; clc
%==========================================================================
% Add paths - edit this section
%==========================================================================

% Project root directory
%projdir = '/gpfs/M2Home/kristina_s/Monash076/Kristina/GenofCog/';
projdir = '/Users/kristinasabaroedin/Documents/UNI/SUBJECTS/PSY4100/';

% Data directory
%datadir = [projdir,'data/'];
datadir = [projdir,'firstlevel/'];

% Analysis output directory
% NOTE, must include .mat file with non imaging data (e.g., GoC_variables.mat)
workdir = [projdir,'secondlevel/'];

% ------------------------------------------------------------------------------
% Load in non imaging data
% ------------------------------------------------------------------------------
cd(workdir)
load('GoC_variables.mat')


% ------------------------------------------------------------------------------
% Remove excluded participants
% ------------------------------------------------------------------------------
B(exclude(:,1),:) = [];
subs(exclude(:,1),:) = [];

clear exclude

numSubs = length(subs);

% ------------------------------------------------------------------------------
% Initialise output directory. If already exists, delete and re-initialise
% ------------------------------------------------------------------------------
WhichROIs = 'TriStri'; % 'TriStri' 'DiMartino'
% Select ROI. 2 = dorsal striatum. 3 = ventral striatum
ROInum = 2;

switch WhichROIs
	case 'TriStri'

		int_dir = 'questionnaires/';

		if ROInum == 1
			ROIname = 'SM'
		elseif ROInum == 2
			ROIname = 'D'
		elseif ROInum == 3
			ROIname = 'V'
		end
	case 'DiMartino'

		int_dir = 'questionnaires_spheres/';

		if ROInum == 1
			ROIname = 'VSs'
		elseif ROInum == 2
			ROIname = 'DC'
		elseif ROInum == 3
			ROIname = 'VRP'
		elseif ROInum == 4
			ROIname = 'DRP'
		elseif ROInum == 5
			ROIname = 'DCP'
		elseif ROInum == 6
			ROIname = 'VSi'
		end
end

% make outdir
outdir = [workdir,int_dir,ROIname,'/'];
if exist(outdir) == 0
	fprintf(1,'Initialising outdir\n')
	mkdir(outdir)
elseif exist(outdir) == 7
	fprintf(1,'Cleaning and re-initialising outdir\n')
	rmdir(outdir,'s')
	mkdir(outdir)
end

conImage = ['con_000',num2str(ROInum),'.img'];

% ------------------------------------------------------------------------------		
% Define covariates
% ------------------------------------------------------------------------------		

Cov_names = {'OLIFE_IA','IQ_4scale','Sex','Age'};
numCov = length(Cov_names);

%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
cd(outdir)

fprintf(1,'Initialising SPM...\n')
spm('defaults','fmri');
spm_jobman('initcfg')

%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.stats.factorial_design.dir = {outdir};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact.name = 'hemi';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact.levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact.dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact.ancova = 0;

% Left hemisphere
	Hemi = 'L';

	switch WhichROIs
		case 'TriStri'
			first_dir = ['/tri/FirstLevel_',Hemi,'/'];
		case 'DiMartino'
			first_dir = ['/spheres/FirstLevel_',Hemi,'_DiMartino/'];
	end

	matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = 1;
	for i = 1:numSubs
		matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans{i} = [datadir,subs{i},first_dir,conImage];
	end

% Right hemisphere
	Hemi = 'R';
	switch WhichROIs
		case 'TriStri'
			first_dir = ['/tri/FirstLevel_',Hemi,'/'];
		case 'DiMartino'
			first_dir = ['/spheres/FirstLevel_',Hemi,'_DiMartino/'];
	end

	matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = 2;
	for i = 1:numSubs
		matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans{i} = [datadir,subs{i},first_dir,conImage];
	end

% Import covariates
	for i = 1:numCov

		matlabbatch{1}.spm.stats.factorial_design.cov(i).cname = Cov_names{i};
		
		matlabbatch{1}.spm.stats.factorial_design.cov(i).iCFI = 1;

		matlabbatch{1}.spm.stats.factorial_design.cov(i).iCC = 1;
		
		% Find index of covariates
		idx = strmatch(Cov_names{i},Bh,'exact');
		
		% duplicate and concatenate to account for hemisphere
		Cov_temp = [B(:,idx);B(:,idx)];

		matlabbatch{1}.spm.stats.factorial_design.cov(i).c = Cov_temp;

	end

matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% ------------------------------------------------------------------------------
% Estimate
% ------------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_est.spmmat = {[outdir,'SPM.mat']};
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% ------------------------------------------------------------------------------
% Contrasts
% ------------------------------------------------------------------------------
matlabbatch{3}.spm.stats.con.spmmat = {[outdir,'SPM.mat']};
matlabbatch{3}.spm.stats.con.delete = 1;

matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'OLIFE_IA+';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [0 0 1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'OLIFE_IA-';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [0 0 -1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';


% ------------------------------------------------------------------------------
% Run
% ------------------------------------------------------------------------------
spm_jobman('run',matlabbatch);
clear matlabbatch

cd(outdir)
