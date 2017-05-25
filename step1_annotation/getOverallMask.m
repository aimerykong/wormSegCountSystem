function imDisplay = getOverallMask(imDisplay, maskDisplay)



TMP = imDisplay(:,:,1); TMP(maskDisplay==1) = 255; imDisplay(:,:,1) = TMP; %255 is white
TMP = imDisplay(:,:,2); TMP(maskDisplay==1) = 0; imDisplay(:,:,2) = TMP;   %0 is black
TMP = imDisplay(:,:,3); TMP(maskDisplay==1) = 0; imDisplay(:,:,3) = TMP;