function Y = GetDesignMatrix(qhat,b)
    images_num = size(qhat,2);
    qx = qhat(1,:);
    qy = qhat(2,:);
    qz = qhat(3,:);
    
    Y = [ones(1,images_num); -b.*qx.^2; -2*b.*qx.*qy; -2*b.*qx.*qz; -b.*qy.^2; -2*b.*qy.*qz; -b.*qz.^2]';
end