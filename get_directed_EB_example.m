
directory = './example_data/';

%make output directory for edge betweenness files
folder = char(strcat(directory, '/EB/'));

if ~exist(folder, 'dir')
   mkdir(folder)
end

%Read network file and convert to Sparse Matrix
f = char(strcat(directory,'10x10grid_edgelist.txt'));
A_sparse = spconvert(...
    table2array(...
        readtable(f)...
    )...
);

%Read node metadata file
nodes = table2array(...
   readtable(...
       char(strcat(directory, '10x10grid_node_metadata.txt'))...
   )...
);

%initialize input vectors
vPop = nodes(:,3);
vStores = nodes(:,4);
vCoords = nodes(:,1:2);
N = size(vStores,1);

[NB, directed_EB] = get_ebw_make_od(A_sparse, N, 1, vStores, vPop, vCoords);

%convert directed to undirected edge betweenness
undirected_EB = triu(directed_EB,1)'+tril(directed_EB);

%scale = 1 so that EB is saved as-is, to be scaled later
scales = [1];

%Save undirected
for scale=scales
    EB_sparse = sparse(undirected_EB*scale);
    [row col v] = find(EB_sparse);
    filename = char(...
        strcat(...
            folder,...
            'eb.txt'...
        )...
    );
    dlmwrite(filename,[col row v], 'delimiter', ' ');
end

%Save upper triangle of directed matrix
for scale=scales
    EB_sparse = sparse(triu(directed_EB*scale,1));
    [row col v] = find(EB_sparse);
    filename = char(...
        strcat(...
            folder,...
            'eb_i.txt'...
        )...
    );
    dlmwrite(filename,[col row v], 'delimiter', ' ');
end

%Save lower triangle of directed matrix
for scale=scales
    EB_sparse = sparse(tril(directed_EB*scale));
    [col row v] = find(EB_sparse);
    filename = char(...
        strcat(...
            folder,...
            'eb_j.txt'...
        )...
    );
    dlmwrite(filename,[col row v], 'delimiter', ' ');
end


