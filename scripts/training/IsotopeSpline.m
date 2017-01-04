% Describe what it does
% Describe the input
% Describe the output
function IsotopeSpline(max_sampled_mass, S, precursor_isotope, infile, outfile_res_hist, outfile_scatter, outfile_res, outfile_gof, outfile_model)
	% For orientation:
	% M(:,1) = probabilities = Y-axis
	% M(:,2) = precursor masses = X-axis
	M = dlmread(infile,'\t',1,0);
	
	
	order = 4;
	
	breakSteps = 250;	% The amount of daltons in between each spline break.
						% We're using the same value for both X and Y.
	min_mass = 50;		% The lower bound of our spline breaks (in daltons)
	
	min_knot = double(idivide(min(M(:,2))-min_mass,int32(breakSteps),'floor')*breakSteps+min_mass); % minimum observed fragment mass rounded
	max_knot = double(idivide(max(M(:,2)),int32(breakSteps),'ceil')*breakSteps+min_mass); % maximum observed fragment mass rounded
	% Create the knots in the proper format
	knots = augknt([min_knot, min_knot:breakSteps:max_knot, max_knot], order);
	% Create spline using least squares approximation in the B-spline format (better for creation)
	sp = spap2(order, M(:,2), M(:,1));
	% Convert from B-spline to piecewise polynomial (pp) format. (better for evaluation)
	sp_pp = fn2fm(sp,'pp');
	sp_pp.breaks
	
	% Create a figure of the scatter plot and spline
	xlabel('precursor mass');
	ylabel('probability');
	title(sprintf(strcat('Precursor isotope: ',precursor_isotope, '\n Sulfurs in peptide: ', S)))
	scatter(M(:,2),M(:,1))
	hold on;
	fnplt(cs, 2, 'g');
	hold off;
	
	[pathstr,name,ext] = fileparts(outfile_scatter);
	print(outfile_scatter,strcat('-d', ext(2:end)));
	
	% Calculate GOF statistics
	[RMSD meanD Rsq residuals] = goodnessOfFitStatistics(M, sp_pp);

	% Write the GOF statistics with model description to a file
	fileID = fopen(outfile_gof,'w');
	fprintf(fileID, '%s', strcat([S, ' ', precursor_isotope, ' ']));
	fprintf(fileID, '%3.5f %3.5f %3.5f\n', [RMSD meanD Rsq]);
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
	print(outfile_res3D,strcat('-d', ext(2:end)));
	
	%writeModelXML(outfile_model, sp_pp, S, precursor_isotope);	
	
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
		residuals(i) = M(i,1) - fnval( sp_pp, M(i,2) ); % M(i,1) is our observed value, and fnval() calculates are predicted value
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
  	if S~='NA'
		fprintf(fileID, strcat([' S=''', S, '''']));
	end
	% Write remaining attributes
	fprintf(fileID, strcat([' PrecursorIsotope=''',precursor_isotope,''' Order=''',num2str(sp_pp.order(1)),'''>\n']));

	% Write the <massBreaks> tag and its attributes
	writeBase64BinaryArrayXML(fileID, 'massBreaks', '32', 'little', num2str(sp_pp.pieces + 1), convertAndEncode(sp_pp.breaks));

	% We initialize the 1-D array that will stores the coefficients to its maximum possible size.
	% The maximum size is the number of patches * the number of coefficients per patch.
	coefs_out = zeros(1, sp_pp.pieces * sp_pp.order );

	% Iterate through the coefficients and store them in the 1-D array.
	ii = 0;
	for i = 1:sp_pp.pieces
        for j = 1:sp_pp.order
			ii = ii+1;
            coefs_out(ii) = sp_pp.coefs(i, j);          
        end
	end

	% Write the <coefficients> tag and its attributes
	writeBase64BinaryArrayXML(fileID, 'coefficients', '32', 'little', num2str(ii), convertAndEncode(coefs_out(1:ii)));

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

	% After testing, double precision was clearly unnecessary.
	% Convert to single precision to half the memory and disk space requirement.
	input = single(input); 
	
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
