% Describe what it does
% Describe the input
% Describe the output
function IsotopeSpline(max_sampled_mass, S, precursor_isotope, infile, outfile_res_hist, outfile_scatter, outfile_res3D, outfile_gof, outfile_model)
	% For orientation:
	% M(:,1) = probabilities = Y-axis
	% M(:,2) = precursor masses = X-axis
	M = dlmread(infile,'\t',1,0);
	
	% Create a figure of the scatter plot
	xlabel('precursor mass');
	ylabel('probability');
	title(sprintf(strcat('Precursor isotope: ',precursor_isotope, '\n Sulfurs in peptide: ', S)))
	scatter(M(:,2),M(:,1))
	[pathstr,name,ext] = fileparts(outfile_scatter);
	print(outfile_scatter,strcat('-d', ext(2:end)));

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
		residuals(i) = M(i,1) - fnval( sp_pp, {M(i,2),M(i,3)} ); % M(i,1) is our observed value, and fnval() calculates are predicted value
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
function writeModelXML(outfile_path, sp_pp, S, CS, Se, CSe, precursor_isotope, fragment_isotope)
	% Open our output file for writing
	fileID = fopen(outfile_path,'w');

	% Open the <model> tag
	fprintf(fileID, '\t<model');
	% Write composition attributes if this was a sulfur specific model
  if S~='NA'
		fprintf(fileID, strcat([' S=''', S, ''' CS=''', CS, ''' Se=''', Se, ''' CSe=''', CSe, '''']));
	end
	% Write remaining attributes
	fprintf(fileID, strcat([' PrecursorIsotope=''',precursor_isotope,''' FragmentIsotope=''',fragment_isotope,''' FragmentOrder=''',num2str(sp_pp.order(1)),''' PrecursorOrder=''',num2str(sp_pp.order(1)),'''>\n']));

	% Write the <fragmentMassBreaks> tag and its attributes
	writeBase64BinaryArrayXML(fileID, 'fragmentMassBreaks', '32', 'little', num2str(length(sp_pp.breaks{1})), convertAndEncode(sp_pp.breaks{1}));

	% Write the <precursorMassBreaks> tag and its attributes
	writeBase64BinaryArrayXML(fileID, 'precursorMassBreaks', '32', 'little', num2str(length(sp_pp.breaks{2})), convertAndEncode(sp_pp.breaks{2}));

	% The format that the coefficients are stored in the tensor product spline is super confusing to me.
	% Reshape them to make them easier to iterate through.
	coefs = reshape(sp_pp.coefs, [1, sp_pp.pieces(1), sp_pp.order(1), sp_pp.pieces(2), sp_pp.order(2)]);
	
	% We want to store the coefficients in a 1-D array before encoding them into Base64
	% We're going to store all the coefficients for a single patch, 
	% followed by the coefficients for next patch, and then the next patch, etc.
	% The coefficients c(1:16) for a particular patch, (assuming the order of the splines is 4 in both x and y)
	% are written in an order such that the equation to evaluate the value for a patch is: 
	% value = c(1)*x^3*y^3 + c(2)*x^3*y^2 + c(3)*x^3*y + c(4)*x^3
	% 		+ c(5)*x^2*y^3 + c(6)*x^2*y^2 + c(7)*x^2*y + c(8)*x^2
	% 		+ c(9)*x*y^3   + c(10)*x*y^2  + c(11)*x*y  + c(12)*x 
	% 		+ c(13)*y^3    + c(14)*y^2    + c(15)*y    + c(16)
	%      
	% Where x = fragment mass - patch's corresponding fragment mass break
	%           
	%      e.g. if a patch is valid for fragment masses = 550-800
	%           and your fragment mass is 650, then x = 650-550 = 100
	%
	%   and y = precursor mass - patch's corresponding precursor mass break
	%
	%      e.g. if a patch is valid for precursor masses = 550-800
	%           and your precursor mass is 650, then y = 650-550 = 100
	%
	%
	%
	% The tensor product spline, sp_pp, contains a full rectangular grid of patches.
	% Our data, however, only covers a triangle because the fragment mass must be smaller than the precursor mass.
	% This means there are many patches (~50%) that are invalid for our model and we won't write them to disk.
	% 
	% We initialize the 1-D array that will stores the coefficients to its maximum possible size.
	% The maximum size is the number of patches * the number of coefficients per patch.
	coefs_out = zeros(1, sp_pp.pieces(1) * sp_pp.pieces(2) * sp_pp.order(1) * sp_pp.order(2));

	% Iterate through the coefficients and store them in the 1-D array.
	ii = 0;
	for i = 1:sp_pp.pieces(2)
	    for j = 1:sp_pp.pieces(1)
			% The patch's fragment mass must be <= precursor mass,
			% otherwise the patch is invalid and we don't write it to disk.
	        if sp_pp.breaks{1}(j) <= sp_pp.breaks{2}(i)
	            for c1 = 1:sp_pp.order(1)
	                for c2 = 1:sp_pp.order(2)
						ii = ii+1;
	                    coefs_out(ii) = coefs(1, j, c1, i, c2);          
	                end
	            end
	        end
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
