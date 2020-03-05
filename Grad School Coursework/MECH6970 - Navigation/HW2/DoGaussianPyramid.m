GaussianPyramid %Run GaussianPyramid first
figure
DoGIs=[];
for i=1:1:pyrLevels-1
    Js = IsPyr(:,(i-1)*size(Is,2)+(1:size(Is,2)))-IsPyr(:,i*size(Is,2)+(1:size(Is,2))); %Computer DoG
    DoGIs = [DoGIs Js];
    subplot(1,pyrLevels-1,i)
    %imshow(histeq(Js)); %Try this!
    imshow(Js,'DisplayRange',[min(min(Js)) max(max(Js))]); %Display DoG ith level
    title(['DoGPyrLevel: ',num2str(i)])
end
