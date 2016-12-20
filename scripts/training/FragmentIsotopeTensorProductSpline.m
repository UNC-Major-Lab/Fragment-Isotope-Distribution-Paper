% Describe what it does
% Describe the input
% Describe the output
function FragmentIsotopeTensorProductSpline(max_sampled_mass,S,CS,Se,CSe,precursor_isotope,fragment_isotope,infile,outfile_spline,outfile_res_hist,outfile_scatter,outfile_res3D,outfile_gof,outfile_model)
	% For orientation:
	% M(:,1) = probabilities = Z-axis
	% M(:,2) = fragment masses = X-axis
	% M(:,3) = precursor masses = Y-axis
	M = dlmread(infile,'\t',1,0);
	
	% Create a figure of the 3D scatter plot
	xlabel('fragment mass');
	ylabel('precursor mass');
	zlabel('probability');
	title(sprintf(strcat('Precursor isotope: ',precursor_isotope, ' Fragment isotope: ', fragment_isotope, '\n Sulfurs in fragment: ', S, ' Sulfurs in complement: ', CS, '\nSeleniums in fragment: ', Se, ' Seleniums in complement: ', CSe)))
	scatter3(M(:,2),M(:,3),M(:,1),1,M(:,1))
	[pathstr,name,ext] = fileparts(outfile_scatter);
	print(outfile_scatter,strcat('-d', ext(2:end)));

	% The input data is in scatter plot form,
	% but it needs to be gridded data to create a tensor product spline.
	%
	% First, determine the bounds of the grid
	breakSteps = 250;	% The amount of daltons in between each spline break.
						% We're using the same value for both X and Y.
	min_mass = 50;		% The lower bound of our spline breaks (in daltons)
	
	min_knot_x = double(idivide(min(M(:,2))-min_mass,int32(breakSteps),'floor')*breakSteps+min_mass); % minimum observed fragment mass rounded
	max_knot_x = double(idivide(max(M(:,2)),int32(breakSteps),'ceil')*breakSteps+min_mass); % maximum observed fragment mass rounded

	min_knot_y = double(idivide(min(M(:,3))-min_mass,int32(breakSteps),'floor')*breakSteps+min_mass); % minimum observed precursor mass rounded
	max_knot_y = double(idivide(max(M(:,3)),int32(breakSteps),'ceil')*breakSteps+min_mass); % maximum observed precursor mass rounded
		
	% These are the x,y coordinates of our grid
	grid_step_size = 10;
	Y = min_knot_y:grid_step_size:max_knot_y;
	X = min_knot_x:grid_step_size:max_knot_x;

	nx = length(X);
	ny = length(Y);

	% Our scatter data does not line up perfectly with the grid,
	% so we need to perform interpolation to determine the values at each grid point
	% We will use the scatteredInterpolant class for this, but it requires the points
	% you want the values of to be in scatter plot format.
	xy_grid_points_in_scatter_format = zeros(nx*ny,2);
	ii = 1;
	for y = Y
		for x = X
			xy_grid_points_in_scatter_format(ii,1) = x;
			xy_grid_points_in_scatter_format(ii,2) = y;
			ii = ii+1;
		end
	end

	% Create the scatteredInterpolant based on our original data
	interpolant = scatteredInterpolant(M(:,2),M(:,3),M(:,1),'linear','nearest');
	% Interpolate the probability values at our grid points (stored in scatter format)
	z_interpolated_in_scatter_format = interpolant(xy_grid_points_in_scatter_format);

	% Now we need to take our interpolated probabilties that are organized for scatter data
	% and convert them to a format for gridded data.
	Z = zeros(nx, ny);

	ii=1;
	for i = 1:length(Y)
		for j = 1:length(X)
			Z(j,i) = z_interpolated_in_scatter_format(ii);
			ii = ii+1;
		end
	end
	
	% Now that we have our input/training data on a grid, 
	% we have to define the parameters of our tensor-product spline
	order = 4;	% The order of the polynomial for each spline piece. 
			  	% We're using the same order for both X and Y dimensions.

	% We will use the same breaks/knots for all the models (different sulfurs, different isotopes).
	% This way, when we're looking at say fragment mass = 500, and precursor mass = 1000, 
	% and we evaluate each of the models with those parameters to see which fits best, then
	% we only need to do one binary search to find the correct index of the patch, because 
	% it'll be the same for all the models.
	%
	% The negative of this is that the breaks aren't at the optimal positions.
	%
	% Not all models will start and end at the same breaks though because for example
	% you can't have a fragment with 6 sulfurs that's only 300Da. So for that model, the
	% the first break would start at a higher number, but it would still be one of the internal
	% breaks of the other models.

	% Create the knots in the proper format
	knots_x = augknt([min_knot_x, min_knot_x:breakSteps:max_knot_x, max_knot_x], order);
	knots_y = augknt([min_knot_y, min_knot_y:breakSteps:max_knot_y, max_knot_y], order);


	sp_y = spap2(knots_y,order,Y,Z);
	coefsy = fnbrk(sp_y,'c');
	sp_x = spap2(knots_x,order,X,coefsy.');
	
	sp_pp_y = fn2fm(sp_y,'pp');
	sp_pp_x = fn2fm(sp_x,'pp');
	
	% Create a tensor product spline in the B-spline format (better for creation)
	sp = spap2({knots_x, knots_y},[order, order], {X,Y}, Z);
	% Convert from B-spline to piecewise polynomial (pp) format. (better for evaluation)
	sp_pp = fn2fm(sp,'pp');

	% To plot a mesh of the spline, we evaluate the spline at the grid break points
	% The values at invalid regions (fragment mass >= precursor mass) must be set
	% to nan so they are ignored during plotting.
	values = fnval( sp_pp, {sp_pp.breaks{1},sp_pp.breaks{2}} );

	for i = 1:length(sp_pp.breaks{1})
	    for j = 1:length(sp_pp.breaks{2})
	        if sp_pp.breaks{1}(i) >= sp_pp.breaks{2}(j)
	            values(i,j) = nan;
	        end
	    end
	end
	
	% Plot mesh of the spline
	mesh(sp_pp.breaks{1},sp_pp.breaks{2},values.');
	xlabel('fragment mass');
	ylabel('precursor mass');
	zlabel('probability');
	title(sprintf(strcat('Precursor isotope: ',precursor_isotope, ' Fragment isotope: ', fragment_isotope, '\n Sulfurs in fragment: ', S, ' Sulfurs in complement: ', CS, '\nSeleniums in fragment: ', Se, ' Seleniums in complement: ', CSe)))                    
	[pathstr,name,ext] = fileparts(outfile_spline);
	print(outfile_spline,strcat('-d', ext(2:end)));

	% Calculate GOF statistics
	[RMSD meanD Rsq residuals] = goodnessOfFitStatistics(M, sp_pp);

	% Write the GOF statistics with model description to a file
	fileID = fopen(outfile_gof,'w');
	fprintf(fileID, '%s', strcat([S,' ',CS,' ',Se,' ',CSe,' ',precursor_isotope,' ',fragment_isotope,' ']));
	fprintf(fileID, '%3.5f %3.5f %3.5f\n', [RMSD meanD Rsq]);
	fclose(fileID);

	% Create histogram of the residuals. Ideally this will be normally distributed and centered at 0.
	figure();
	hist(residuals,50);
	xlabel('residual');
	title(sprintf(strcat('Precursor isotope: ',precursor_isotope, ' Fragment isotope: ', fragment_isotope, '\n Sulfurs in fragment: ', S, ' Sulfurs in complement: ', CS, '\nSeleniums in fragment: ', Se, ' Seleniums in complement: ', CSe)))
	[pathstr,name,ext] = fileparts(outfile_res_hist);
	print(outfile_res_hist,strcat('-d', ext(2:end)));
	
	% Create 3D scatterplot of residuals. This will show us any hot spots of high error.
	% The PDF takes up a lot of space but looks better.
	% If this isn't for publication, then output as a non-vector graphics format
	figure();
	scatter3(M(:,2),M(:,3),residuals,1,residuals)
	xlabel('fragment mass');
	ylabel('precursor mass');
	zlabel('residual');
	title(sprintf(strcat('Precursor isotope: ',precursor_isotope, ' Fragment isotope: ', fragment_isotope, '\n Sulfurs in fragment: ', S, ' Sulfurs in complement: ', CS, '\nSeleniums in fragment: ', Se, ' Seleniums in complement: ', CSe)))
	[pathstr,name,ext] = fileparts(outfile_res3D);
	print(outfile_res3D,strcat('-d', ext(2:end)));

	writeModelXML(outfile_model,sp_pp,S,CS,Se,CSe,precursor_isotope,fragment_isotope);	

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
		fprintf(fileID, strcat([' S=''',S,''' CS=''',CS,''' Se=''',Se,''' CSe=''',CSe,''']);
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
