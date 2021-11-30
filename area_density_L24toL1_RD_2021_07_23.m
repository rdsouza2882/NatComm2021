%%% get density of L1 and L24 across multiple pathways
% updated 2020/12/13 on thermaltake

%% set parameters
filelist_excel = 'C:\Users\Burkhalter Lab\Documents\Rinaldo\Optical Density\2020_analyses_file_list_2021-07-23.xlsx'; 
%filelist_excel = 'C:\Users\Burkhalter Lab\Documents\Rinaldo\Optical Density\2020_analyses_file_list_2020-12-08_forFig1.xlsx';
%filelist_excel = 'C:\Users\Burkhalter Lab\Documents\Rinaldo\Optical Density\2020_analyses_file_list_V1_AL_forpixelhistograms.xlsx'; 
   excel_rows_to_analyze = []; % if empty, analyze all rows
%  excel_rows_to_analyze = [2:15];
group_varnames(1).name = 'L1';
group_varnames(2).name = 'L24';
black_white_file_extension = '.jpg';
proj_file_extension = '.tif';
roi_directory_varname = 'roi_directory';
proj_directory_varname = 'file_directory'; 
output_analyzed_pixels = 1; % if true, output all (pre-baselined) analyzed pixels for each case in the filelist table
areaDensity_ops.area_proportion = 1; % takes only the top specified fraction of pixels from the ROI to analyze(set = 1 to not use this option)
areaDensity_ops.intensity_vals_proportion = 0.70; % take only the top specified fraction of pixel intensity vals to analyze (set = 1 to not use this option) (percentage of max pixel value)
areaDensity_ops.do_baseline_subtraction = 1; % if true, include the step of baseline-subtracting intensity values
areaDensity_ops.separate_projfile_directory = true; % toggle wheter proj files are in different directory from roi and baseline files
areaDensity_ops.save_analyzed_image = 0; % if true, save the analyzed image in the output table (may be useful for looking at blurred images)
areaDensity_ops.blur_image = 1; % if true, blur image before analyzing pixels
    areaDensity_ops.diskblurradius_um = 1; % blur radius; the value used used in findpatches.m = 29um

% loadbw options
areaDensity_ops.loadbw_ops.allow_multiple_nonzero = true;
% areaDensity_ops.loadbw_ops.logical_zero_value = [16]; % in this dataset, 16 = ROI
areaDensity_ops.loadbw_ops.reverse_black_white = 0; % flip zeros and nonzeros




%% run analysis
filelist = readtable(filelist_excel);
if ~isempty(excel_rows_to_analyze)
    filelist = filelist(excel_rows_to_analyze-1,:); % keep only specified rows
end
filelist.Properties.VariableNames{strcmp(filelist.Properties.VariableNames, 'xCase')} = 'case'; % correct irregular varname
filelist = movevars(filelist, 'case', 'After', 'pathway');
filelist = filelist( ~cellfun(@(x)isempty(x),filelist.proj_file) ,:); % delete excel rows without a listed proj file

pathway_list = unique(filelist.pathway, 'stable'); %for grouping individual pathways (2 cases per pathway)
cases_list = unique(filelist.uniquepathway, 'stable'); %for grouping individual cases (2-4 sections per case)

npathways = length(pathway_list); 
ncases = length(cases_list);

n_sorting_groups = length(group_varnames);
nfiles = height(filelist); 
areaDensity_ops.output_analyzed_pixels = output_analyzed_pixels; % set option to output or not output analyzed pixels

% format proj file names
if areaDensity_ops.separate_projfile_directory % add directory if necessary
    filelist.proj_file = strcat(filelist{:,proj_directory_varname}, filesep, filelist.proj_file);
end
filelist.proj_file = strcat(filelist.proj_file, proj_file_extension); % add extension    

for igroup = 1:n_sorting_groups
    thisgroup = group_varnames(igroup).name;
    group_varnames(igroup).roi = [thisgroup, '_file'];
    group_varnames(igroup).baseline = [thisgroup, 'baseline_file'];
    filelist_temp = filelist;

    filelist_temp.directory = filelist{:,roi_directory_varname};
    filelist_temp.exposure_sec = NaN(nfiles,1);
    % specify group to analyze, format input
    filelist_temp.roi_file = strcat(filesep, filelist{:,group_varnames(igroup).roi}, black_white_file_extension); 
    filelist_temp.baseline_file = strcat(filesep, filelist{:,group_varnames(igroup).baseline}, black_white_file_extension); 
    filelist_temp = areaDensity(filelist_temp, areaDensity_ops); % run areaDensity on this group
    filelist{:,['intensPerPix_', thisgroup]} = filelist_temp.intensPerPix; % move results to main table
    if areaDensity_ops.do_baseline_subtraction; filelist{:,['baseline_', thisgroup]} = filelist_temp.baseline; end % move results to main table; 
    if areaDensity_ops.save_analyzed_image; filelist{:,['analyzed_image_', thisgroup]} = filelist_temp.analyzed_image; end % move results to main table;
    if output_analyzed_pixels
        filelist{:,['analyzed_pix_', thisgroup]} = filelist_temp.analyzed_pix; % move results to main table
        if areaDensity_ops.do_baseline_subtraction
            for ifile = 1:height(filelist_temp)
                filelist{ifile,['analyzed_pix_baselined_', thisgroup]} = {filelist_temp.analyzed_pix{ifile} - filelist_temp.baseline(ifile)}; % subtract baselines
            end
        end
    end
    filelist = movevars(filelist, {['intensPerPix_', thisgroup]}, 'After', 'case');
    if areaDensity_ops.do_baseline_subtraction; filelist = movevars(filelist, {['baseline_', thisgroup]}, 'After', ['intensPerPix_', thisgroup]); end
   clear filelist_temp 
end

% get ratio between pixel intensities in group 1 and group 2
ratio_varname1 = ['ratio_', group_varnames(2).name, '_', group_varnames(1).name];
ratio_varname2 = ['ratio_', group_varnames(2).name, '_', group_varnames(1).name, 'plus', group_varnames(2).name];
%ratio_varname3 = ['ratio_', group_varnames(1).name, '_', group_varnames(2).name];
%ratio_varname4 = ['ratio_', group_varnames(1).name, '_', group_varnames(1).name, 'plus', group_varnames(2).name];
filelist{:,ratio_varname1} = filelist{:,['intensPerPix_', group_varnames(2).name]} ./ filelist{:,['intensPerPix_', group_varnames(1).name]};
filelist{:,ratio_varname2} = filelist{:,['intensPerPix_', group_varnames(2).name]} ./ (filelist{:,['intensPerPix_', group_varnames(1).name]} + filelist{:,['intensPerPix_', group_varnames(2).name]});
%filelist{:,ratio_varname3} = filelist{:,['intensPerPix_', group_varnames(1).name]} ./ filelist{:,['intensPerPix_', group_varnames(2).name]};
%filelist{:,ratio_varname4} = filelist{:,['intensPerPix_', group_varnames(1).name]} ./ (filelist{:,['intensPerPix_', group_varnames(1).name]} + filelist{:,['intensPerPix_', group_varnames(2).name]});
filelist = movevars(filelist, ratio_varname2, 'After', 'case');

%% group results by pathway
nans = NaN(size(pathway_list));
group_results = table(pathway_list, cell(size(pathway_list)), nans, nans, nans, nans, 'VariableNames', {'pathway',ratio_varname2,'nfiles','mean','median','sem'});
for ipath = 1:length(pathway_list)
    thesevals = filelist{strcmp(filelist.pathway,pathway_list{ipath}), ratio_varname2};
    group_results{ipath,ratio_varname2} = {thesevals};
    group_results.nfiles(ipath) = length(thesevals);
    group_results.mean(ipath) = mean(thesevals);
    group_results.sem(ipath) = std(thesevals) ./ sqrt(length(thesevals));
    group_results.median(ipath) = nanmedian(thesevals);
end

%% group results by case number
nans = NaN(size(cases_list));
group_results = table(cases_list, cell(size(cases_list)), nans, nans, nans, nans, 'VariableNames', {'case number',ratio_varname2,'nfiles','mean','median','sem'});
for icase = 1:length(cases_list)
    thesevals = filelist{strcmp(filelist.uniquepathway, cases_list{icase}), ratio_varname2};
    group_results{icase,ratio_varname2} = {thesevals};
    group_results.nfiles(icase) = length(thesevals);
    group_results.mean(icase) = mean(thesevals);
    group_results.sem(icase) = std(thesevals) ./ sqrt(length(thesevals));
    group_results.median(icase) = nanmedian(thesevals);
end

group_results


% for ipath = 1:npathways
% thispath = pathway_list{ipath};
% matchrows = strcmp(filelist.pathway, thispath); 
    
    