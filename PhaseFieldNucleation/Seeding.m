function [phi]=Seeding(phi,C3SThick,Dcsh,Point,index)

phi(C3SThick:C3SBound,Point(index)-Dcsh:Point(index)+Dcsh,2)=0.8;  

end

