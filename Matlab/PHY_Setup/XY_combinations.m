vec1 = [-24 -21 -18 -15 -12 -9 -6 -3 0 3];
vec2 = [-12 -9 -6 -3 0 3 6 9 12];

[p,q] = meshgrid(vec1, vec2);
pairs = [p(:) q(:)];


comb = [];
for j = 1:100

  comb = [comb; pairs(randperm(length(pairs)),:)];

end

pth             = '/Users/fschneider/Desktop/';
fid             = fopen([pth 'x_position.txt'],'w');
vec             = sprintf('%f\n', comb(:,1));                 % Convert to tab-spaced vector
fprintf(fid, vec);                                                      % Write to file
fclose(fid);                                                            % Close file

pth             = '/Users/fschneider/Desktop/';
fid             = fopen([pth 'y_position.txt'],'w');
vec             = sprintf('%f\n', comb(:,2));                 % Convert to tab-spaced vector
fprintf(fid, vec);                                                      % Write to file
fclose(fid);                                                               % Close file
