clear all;
close all;

x = -3:3;
y = fliplr(-3:3);
[x_mat, y_mat]  = meshgrid( x,y);

z = x_mat + 1i*y_mat;

z_vect = z(:);

dist_z_to_z = abs( repmat(z_vect,1,49 ) - repmat(z_vect.',49,1 )  );

PL_dB = 10*log10(dist_z_to_z.^( -3 ));


        