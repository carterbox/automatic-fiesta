function data = calc_mean_strains(data)
% Calculates mean strains from the strains

Exx = data.getResponse('Exx');
Eyy = data.getResponse('Eyy');
Ezz = data.getResponse('Ezz');
Exy = data.getResponse('Exy');
Exz = data.getResponse('Exz');
Eyz = data.getResponse('Eyz');

[Em,~,~,~] = mean_strains(Exx,Eyy,Ezz,Eyz,Exy,Exz);

data = data.addResponse('Em', Em, ''); 
end

    