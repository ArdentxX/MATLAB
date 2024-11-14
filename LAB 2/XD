close all;

x= 52;
y= 38;
kwadratx = 26;
kwadraty = 17;
rozmx =1;
rozmy =1;

nadajniki = [12.5 26.5 40.5;2.5 2.5 2.5];
odbiorniki =[8 30 44.5;36.5 36.5 36.5];

plot(x,y)
hold on
% blocked =zeros()
for i =1:3
    for j =1:3
    k = wektorsektor(nadajniki(1,i),nadajniki(2,i),odbiorniki(1,j),odbiorniki(2,j),kwadratx,kwadraty,rozmx,rozmy);
    if k == -1
        line([nadajniki(1,i),odbiorniki(1,j)],[nadajniki(2,i),odbiorniki(2,j)],"Color","g","LineWidth",4)
    end
    if k == 1
        line([nadajniki(1,i),odbiorniki(1,j)],[nadajniki(2,i),odbiorniki(2,j)],"Color","r","LineWidth",4)

    end

    end
end

for i =0:x
    for j=0:y
        rectangle('Position',[i j 1 1])

    end
end



plot(nadajniki(1,:),nadajniki(2,:),"r*")
plot(odbiorniki(1,:),odbiorniki(2,:),"b*")
rectangle('Position',[kwadratx kwadraty rozmx rozmy],'FaceColor',"b")
ylim([0,38])
xlim([0,52])
grid on
hold off
