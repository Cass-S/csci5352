function [max_z, max_q, outinfo] = community_modularity(A, varargin)

% function [z, q, outinfo] = community_modularity(A, optional_inputs)
%
% feature: This function conduct the greedy agglomerative algorithm to find
%          community structure that maximizes the network's modularity (Q).
%
% input:   A      adjacency matrix
%
% output:  max_z    a vector of community membership that maximize Q
%          max_q    the maximum modularity value
%          outinfo  saves max_z,max_q, and z and q for each merging step
%
% options: 'vis', 'visualize'   visualizes the network using the webweb
%                               tool (http://danlarremore.com/webweb/)
%
% example:
%       A = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; ...
%            0 0 0 1 0 1; 0 0 0 1 1 0];
%       [z, q, outinfo] = greedy_agglo(A)  % returns max z and q, and info
%       [z, q, outinfo] = greedy_agglo(A, 'vis')
%
% karate: A = [0 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 0 ;1 0 1 1 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 1 0 0 0 ;1 1 0 1 0 0 0 1 1 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 0 ;1 1 1 0 0 0 0 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1 % 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 ;0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ;1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ;0 0 0 0 0 0 0 0 % 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 ;0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 ;1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 ;1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 % 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 1 1 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 1 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 ;0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 ;0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 % 0 0 0 0 0 0 0 0 0 1 0 1 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1 1 ;0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 ;1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 0 0 1 1 ;0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 1 0 0 1 0 1 0 1 1 0 0 0 0 0 1 1 1 0 1 ;0 0 0 0 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 0 0 1 0 0 1 1 1 1 1 1 1 0 ]
%
%
% All calculations are based on the lecture note of Aaron Clauset's Network
% analysis and modeling class (Fall 2014).
% see  http://tuvalu.santafe.edu/~aaronc/courses/5352/


dovis = false;
dogiven = false;

% Parse the optional inputs
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'vis', 'visualize'} % visualize the network (using Webweb)
                dovis = true;
                originalA = A; % keep the original A for its visualization
            case {'given_solution', 'given'} % visualize the known category
                dogiven = true;
                z_given = varargin{i+1};
        end
    end
end

m = sum(sum(A))./2;

%% 1. initialize z, q values
z = 1:size(A,1);                    % Each vertex is in its own group.
q = modularity_wani(A, z);          % Calculate the initial Q value

%% delta_q function handle
delta_q = @(A,m) A./m - 2.*(sum(A)'*sum(A)) .* (A~=0) ./ ((2*m)^2);

                    % Eq.(8) of Lecture note 5. In this implementation,
                    % A./m is 2*e_ij in Eq.(8), and ki./(2*m) is a_i.
                    % Multiplying (A~=0) is added because this calculation
                    % should be conducted only on the existing pairs.
                    % This value has been validated using a test against
                    % diff(q). They were same.

%% 2. greedy agglomeration

i = 1;
while size(A,1) > 2

    dq = delta_q(A,m);
    dq(eye(size(dq))~=0)=-Inf;

    % get the indices for the pair that has the maximum delta q
    [~, idx] = max(dq);
    [~, idx2] = max(max(dq));
    q_delta(i) = max(max(dq));

    % merge the network: see the subfunction "merge"
    dimensions = size(z);
    if dimensions(1) <= 1
        [A, z(i+1,:)] = merge(A, z, sort([idx(idx2) idx2]));
    else
        [A, z(i+1,:)] = merge(A, z(i,:), sort([idx(idx2) idx2]));
    end
    q(i+1) = q(i) + q_delta(i);
                    % calculate the modularity of the merged network.
    i = i + 1;
end

[max_q, maxidx] = max(q);
max_z = z(maxidx,:);

% save output values
outinfo.max_z = max_z;
outinfo.max_q = max_q;
outinfo.z = z;
outinfo.q = q;
outinfo.q_delta = q_delta;

%% 3. visualization using webweb (optional)
if dovis
    webwebdir = '/Users/cassi/Documents/MATLAB/webweb';
    try
        cd(webwebdir);
    catch
        disp('Please specify the directory that includes the Webweb tool in the code.');
        return
    end
    ww.display.w = 300;
    ww.display.h = 300;
    % Increase the charge and the gravity
    ww.display.c = 40;
    ww.display.g = 0.05;
    ww.display.l = 13;
    ww.display.r = 4;
    % Give the file a name
    ww.display.name = 'community_structure_modularity';
    % Name the nodes

    for i=1:size(originalA,1)
        ww.display.metadata.nodeNames{i} = num2str(i);
    end
    % Convert from the matrix to an edgeList.
    [from,to,weight] = find(triu(originalA));
    % Place this edgeList into a network within the webweb struct
    ww.networks.network.edgeList = [from,to,weight];
    alphabeticallity_keys = {'1','2','3','4'};
    ww.display.metadata.alphabeticallity.values = max_z;
    ww.display.metadata.alphabeticallity.categories = alphabeticallity_keys;
    ww.display.colorBy = 'alphabeticallity';

    webweb(ww);
    
    webwebdir = '/Users/cassi/Dropbox/CSCI5352/project';
    cd(webwebdir);

end

end

% ------------------------- SUBFUNCTION -------------------------
function [mA, mz] = merge(A,z,ij)

% This subfunction is merging the vertex pairs that has the largest
% q_delta value.

i = ij(1);
j = ij(2);

mA = A;

mA(i,:) = mA(i,:) + mA(j,:); mA(j,:) = [];
mA(:,i) = mA(:,i) + mA(:,j); mA(:,j) = [];

z2 = unique(z);
ii = z2(i);
jj = z2(j);

mz = z;
mz(z == jj) = mz(ii);

end
