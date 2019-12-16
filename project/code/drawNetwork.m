function q = drawNetwork(A,z,varargin)

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

    for i=1:size(A,1)
        ww.display.metadata.nodeNames{i} = num2str(i);
    end
    % Convert from the matrix to an edgeList.
    [from,to,weight] = find(triu(A));
    % Place this edgeList into a network within the webweb struct
    ww.networks.network.edgeList = [from,to,weight];
    alphabeticallity_keys = {'1','2'};
    ww.display.metadata.alphabeticallity.values = z;
    %ww.display.metadata.alphabeticallity.categories = alphabeticallity_keys;
    ww.display.colorBy = 'alphabeticallity';

    webweb(ww);

    webwebdir = '/Users/cassi/Dropbox/CSCI5352/project';
    cd(webwebdir);
    
end