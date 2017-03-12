function data = calc_von_mises(data)
% Calculate Von Mises strain

Exx = data.getResponse('Exx');
Eyy = data.getResponse('Eyy');
Ezz = data.getResponse('Ezz');
Exy = data.getResponse('Exy');
Exz = data.getResponse('Exz');
Eyz = data.getResponse('Eyz');

Evm = von_mises(Exx,Eyy,Ezz,Eyz,Exy,Exz);

data = data.addResponse('Evm', Evm, ''); 
end

function Evm = von_mises(Exx,Eyy,Ezz,Eyz,Exy,Exz)
% Calculate Von Mises strain

Evm = zeros(size(Exx));

parfor i = 1:numel(Exx)
    
    S = [Exx(i),Exy(i),Exz(i);Exy(i),Eyy(i),Eyz(i);Exz(i),Eyz(i),Ezz(i)]; 
    
    Evm(i) = sqrt(3/2 * trace(S.^2));

end

end