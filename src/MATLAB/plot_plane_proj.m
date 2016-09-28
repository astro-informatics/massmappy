function [] = plot_plane_proj( f, L, a, c)
%UNTITLED Summary of this function goes here

dec = 90 - (a * 180 / L );
ra = c * 180 /(2*L-1);


b = imagesc(ra, dec, real(f(a,c)));
set(gca, 'YDir', 'normal');
set(b,'AlphaData',~isnan(f(a,c)))
xlabel('RA')
ylabel('Dec')
grid

end

