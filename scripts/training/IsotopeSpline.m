% Describe what it does
% Describe the input
% Describe the output
function IsotopeSpline(breakSteps, S, precursor_isotope, infile, outfile_res_hist, outfile_scatter, outfile_res, outfile_gof, outfile_model, outfile_spline_eval)
	% For orientation:
	% M(:,1) = probabilities = Y-axis
	% M(:,2) = precursor masses = X-axis
	M = dlmread(infile,'\t',1,0);
	
	
	order = 4;

	min_mass = 50;		% The lower bound of our spline breaks (in daltons)
	
	min_knot = min(M(:,2));
	max_knot = max(M(:,2));
	% Create uniformly spaced knots in the proper format
	knots = augknt([min_knot, min_knot:breakSteps:max_knot, max_knot], order, 1);
	% Create spline using least squares approximation in the B-spline format (better for creation)
	sp = spap2(knots, order, M(:,2), M(:,1));
	% Optimize knot selection
	sp = spap2(newknt(sp), order, M(:,2), M(:,1));	
	% Convert from B-spline to piecewise polynomial (pp) format. (better for evaluation)
	sp_pp = fn2fm(sp,'pp');
	sp_pp.breaks
	num_neg = testNonNegative(sp_pp, min_knot, max_knot)
	
	% Create a figure of the scatter plot and spline
	xlabel('precursor mass');
	ylabel('probability');
	title(sprintf(strcat('Precursor isotope: ',precursor_isotope, '\n Sulfurs in peptide: ', S)))
	scatter(M(:,2),M(:,1))
	hold on;
	fnplt(sp, 2, 'g');
	hold off;
	
	[pathstr,name,ext] = fileparts(outfile_scatter);
	print(outfile_scatter,strcat('-d', ext(2:end)));
	
	% Calculate GOF statistics
	[RMSD meanD Rsq residuals] = goodnessOfFitStatistics(M, sp_pp);

	% Write the GOF statistics with model description to a file
	fileID = fopen(outfile_gof,'w');
	fprintf(fileID, '%s', strcat([S, ' ', precursor_isotope, ' ']));
	fprintf(fileID, ' %3.5f %3.5f %3.5f %d\n', [RMSD meanD Rsq num_neg]);
	fclose(fileID);
	
	% Create histogram of the residuals. Ideally this will be normally distributed and centered at 0.
	figure();
	hist(residuals, 50);
	xlabel('residual');
	title(sprintf(strcat('Precursor isotope: ',precursor_isotope, '\n Sulfurs in fragment: ', S)))
	[pathstr,name,ext] = fileparts(outfile_res_hist);
	print(outfile_res_hist,strcat('-d', ext(2:end)));
	
	% Create a scatterplot of residuals. This will show us any hot spots of high error.
	% The PDF takes up a lot of space but looks better.
	% If this isn't for publication, then output as a non-vector graphics format
	figure();
	scatter(M(:,2), residuals, 1, residuals)
	xlabel('precursor mass');
	ylabel('residual');
	title(sprintf(strcat('Precursor isotope: ',precursor_isotope,' Sulfurs in fragment: ', S)))
	[pathstr,name,ext] = fileparts(outfile_res);
	print(outfile_res,strcat('-d', ext(2:end)));
	
	writeModelXML(outfile_model, sp_pp, S, precursor_isotope);	
	
	evaluateAndOutputSpline(outfile_spline_eval, sp_pp, min_knot, max_knot);
	
	exit;
end

% This function calculates the residuals, R^2, Root Mean Squared Deviation (RMSD), and average absolute deviation
% for our spline model compared to the data it was trained against. We're not worried about overfitting
% because there is no error associated with our measurements. They should be a good representation of
% the average values we'll see in the real world, but we test that against tryptic peptides in a different program.
function [RMSD meanD Rsq residuals] = goodnessOfFitStatistics(M, sp_pp)
	
	% Calculate residual for each data point.
	residuals = zeros(1,length(M(:,1)));
	for i = 1:length(M(:,1))
		residuals(i) = M(i,1) - fnval( sp_pp, M(i,2) ); % M(i,1) is our observed value, and fnval() calculates the predicted value
	end

	% Calculate the GOF statisitics
	meanD = mean(abs(residuals));
	RMSD = sqrt(mean(residuals.^2));
	average = mean(M(:,1));
	SStot = sum((M(:,1)-average).^2);
	SSres = sum(residuals.^2);
	Rsq = 1 - (SSres/SStot);	
end

% Takes a tensor-product spline model in pp form and writes a partial XML file to disk.
% This will write the information for a single <model> tag.
function writeModelXML(outfile_path, sp_pp, S, precursor_isotope)
	% Open our output file for writing
	fileID = fopen(outfile_path,'w');

	% Open the <model> tag
	fprintf(fileID, '\t<model');
	% Write composition attributes if this was a sulfur specific model
  	if ~strcmp(S,'-1')
		fprintf(fileID, strcat([' S=''', S, '''']));
	end
	% Write remaining attributes
	fprintf(fileID, strcat([' isotope=''',precursor_isotope,''' order=''',num2str(sp_pp.order(1)),'''>\n']));

	% Write the <massBreaks> tag and its attributes
	writeBase64BinaryArrayXML(fileID, 'knots', '64', 'little', num2str(sp_pp.pieces + 1), convertAndEncode(sp_pp.breaks));

	% We initialize the 1-D array that will stores the coefficients to its maximum possible size.
	% The maximum size is the number of patches * the number of coefficients per patch.
	coefs_out = zeros(1, sp_pp.pieces * sp_pp.order );

	% Iterate through the coefficients and store them in the 1-D array.
	ii = 0;
	for i = 1:sp_pp.pieces
		for j = sp_pp.order:-1:1
			ii = ii+1;
			coefs_out(ii) = sp_pp.coefs(i, j);          
		end
	end

	% Write the <coefficients> tag and its attributes
	writeBase64BinaryArrayXML(fileID, 'coefficients', '64', 'little', num2str(ii), convertAndEncode(coefs_out(1:ii)));

	% Close the <model> tag
	fprintf(fileID, '\t</model>\n'); 
	
	% Close the file
	fclose(fileID);
end


% Writes the base64BinaryArrayXML tag and its attributes
function writeBase64BinaryArrayXML(fileID, tag, precision, endian, length, data)

	% Open the tag and write its attributes
	fprintf(fileID, strcat(['\t\t<', tag, ' precision=''', precision, ''' endian=''', endian, ''' length=''', length, '''>']));
	% Write the actual values
	fprintf(fileID, data);
	% Close the tag
	fprintf(fileID, strcat(['</', tag, '>\n']));
end



function output = convertAndEncode(input)
	[str,maxsize,endian] = computer;
	
	% Executive decision to store the data as little-endian
	if endian == 'R'
		input = swapbytes(input); % I actually haven't tested this to be sure it works, but it should...
	end
	
	output = typecast(input,'uint8'); % Convert to byte array
	output = base64encode(output);
	
	% The version of Java's Base64 codec used with Matlab always adds CRLFs every 80 characters or so.
	% The Xerces-C XML parser doesn't get rid of the CRLFs when parsing Base64 data.
	% This messed up decoding of the Base64 data, so get rid of the CRLFs before writing to disk.
	output = regexprep(output, '\r\n', ''); 
end

function output = base64encode(input)
%BASE64ENCODE Encode a byte array using Base64 codec.
%
%    output = base64encode(input)
%
% The function takes a char, int8, or uint8 array INPUT and returns Base64
% encoded string OUTPUT. JAVA must be running to use this function. Note
% that encoding doesn't preserve input dimensions.
%
% See also base64decode

	error(nargchk(1, 1, nargin));
	error(javachk('jvm'));
	if ischar(input)
		input = uint8(input); 
	end

	output = char(org.apache.commons.codec.binary.Base64.encodeBase64Chunked(input))';
end


function output = testNonNegative(sp_pp, min_mass, max_mass)
	output = 0; 
	for i = min_mass:max_mass
		y = fnval( sp_pp, i );
		if y < 0
			output = output+1;
		end
	end
end

function evaluateAndOutputSpline(outfile_path, sp_pp, min_mass, max_mass)
	% Open our output file for writing
	fileID = fopen(outfile_path,'w');
	
	fprintf(fileID, 'precursor.mass\tprobability');
	
	for i = min_mass:max_mass
		y = fnval( sp_pp, i );
		fprintf(fileID, '%d\t%d', i, y);
	end
end

