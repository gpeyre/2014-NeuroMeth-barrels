function Itr = translate_image(I,tx,ty)
    [N,M] = size(I);

    if tx>=0
        idx = [M-tx:M,1:M-tx-1];
    else
        idx = [1-tx:M,1:tx];
    end
    if ty>=0
        idy = [N-ty:N,1:N-ty-1];
    else
        idy = [1-ty:N,1:ty];
    end

    Itr = I(idy,idx);
end