function [gLSSS_N,gLSSS_Np1,gLm1SSS_N,gLm1SSS_Np1] = vrr_1(a,b,c,d,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,Lmax,order,Boys_Table,nz)

if L1 == 1 %PSSS
    
    if order == 0
    g_SSSS_N = g_SSSS_N_2(Kab,Kcd,p,q,alpha,RPQ2,2,Boys_Table); %contains orders 0,1,2
    else
    g_SSSS_N = g_SSSS_N_2(Kab,Kcd,p,q,alpha,RPQ2,2+order,Boys_Table);
    g_SSSS_N = g_SSSS_N(order+1:2+order+1);
    end
    
    %gLSSS_N = RPA*g_SSSS_N(1)-a/p*RPQ*g_SSSS_N(2);
    gLSSS_N = RPA*g_SSSS_N(1)+RWP*g_SSSS_N(2);
    %april 7th 2017. I want to see if defining empty matrices is slower
    %than setting the output to zero
    gLSSS_Np1 = RPA*g_SSSS_N(2)+RWP*g_SSSS_N(3);
    gLm1SSS_N = g_SSSS_N(1);
    gLm1SSS_Np1 = g_SSSS_N(2);
    
elseif L1 == 2 %DSSS

    if order == 0
    g_SSSS_N = g_SSSS_N_2(Kab,Kcd,p,q,alpha,RPQ2,3,Boys_Table); %contains orders 0,1,2,3
    else
    g_SSSS_N = g_SSSS_N_2(Kab,Kcd,p,q,alpha,RPQ2,3+order,Boys_Table);
    g_SSSS_N = g_SSSS_N(order+1:3+order+1);
    end
    
    %g_SSSS_N = g_SSSS_N_2(Kab,Kcd,p,q,alpha,RPQ2,2,Boys_Table); %contains orders 0,1,2
    
    gPSSS_0 = RPA*g_SSSS_N(1)+RWP*g_SSSS_N(2);
    gPSSS_1 = RPA*g_SSSS_N(2)+RWP*g_SSSS_N(3);
    gPSSS_2 = RPA*g_SSSS_N(3)+RWP*g_SSSS_N(4);
    
    %Alternative on the expand function, where the order of expansion is
    %changed.
%     gLSSS_N = [RPA(1)*gPSSS_0;RPA(2)*gPSSS_0(end-L1+1:end);RPA(3)*gPSSS_0(end)]...
%                 +[RWP(1)*gPSSS_1;RWP(2)*gPSSS_1(end-L1+1:end);RWP(3)*gPSSS_1(end)];
%      gLSSS_N = [RPA*gPSSS_0(1);RPA(end-L1+1:end)*gPSSS_0(2);RPA(end)*gPSSS_0(end)]...
%                  +[RWP*gPSSS_1(1);RWP(end-L1+1:end)*gPSSS_1(2);RWP(end)*gPSSS_1(end)];
                       
%     gLSSS_Np1 = [RPA(1)*gPSSS_1;RPA(2)*gPSSS_1(end-L1+1:end);RPA(3)*gPSSS_1(end)]...
%                 +[RWP(1)*gPSSS_2;RWP(2)*gPSSS_2(end-L1+1:end);RWP(3)*gPSSS_2(end)];        
%     gLSSS_Np1 = [RPA*gPSSS_1(1);RPA(end-L1+1:end)*gPSSS_1(2);RPA(end)*gPSSS_1(end)]...
%                  +[RWP*gPSSS_2(1);RWP(end-L1+1:end)*gPSSS_2(2);RWP(end)*gPSSS_2(end)];
            %This is for DSSS functions
     Lm2Terms = 1/2/p*(g_SSSS_N(1)-q/(p+q)*g_SSSS_N(2));        

    
     gLSSS_N = [RPA(1)*gPSSS_0(1)+RWP(1)*gPSSS_1(1)+Lm2Terms
                RPA(2)*gPSSS_0(1)+RWP(2)*gPSSS_1(1)
                RPA(3)*gPSSS_0(1)+RWP(3)*gPSSS_1(1)
                RPA(2)*gPSSS_0(2)+RWP(2)*gPSSS_1(2)+Lm2Terms
                RPA(3)*gPSSS_0(2)+RWP(3)*gPSSS_1(2)
                RPA(3)*gPSSS_0(3)+RWP(3)*gPSSS_1(3)+Lm2Terms];
    
    Lm2Terms_p1 = 1/2/p*(g_SSSS_N(2)-q/(p+q)*g_SSSS_N(3));        

   gLSSS_Np1 = [RPA(1)*gPSSS_1(1)+RWP(1)*gPSSS_2(1)+Lm2Terms_p1
                RPA(2)*gPSSS_1(1)+RWP(2)*gPSSS_2(1)
                RPA(3)*gPSSS_1(1)+RWP(3)*gPSSS_2(1)
                RPA(2)*gPSSS_1(2)+RWP(2)*gPSSS_2(2)+Lm2Terms_p1
                RPA(3)*gPSSS_1(2)+RWP(3)*gPSSS_2(2)
                RPA(3)*gPSSS_1(3)+RWP(3)*gPSSS_2(3)+Lm2Terms_p1];
    
    gLm1SSS_N = gPSSS_0;
    gLm1SSS_Np1 = gPSSS_1;
    
elseif L1 == 3 %FSSS
    %This is the last of the handwritten expressions. So
    %for L1 > 3, I will need to call the function recursively. In that case
    %the function output may need to be of order 0 or higher, and I need to
    %implement an if condition.
    if order == 0
    g_SSSS_N = g_SSSS_N_2(Kab,Kcd,p,q,alpha,RPQ2,4,Boys_Table); %contains orders 0,1,2,3,4
    else
    g_SSSS_N = g_SSSS_N_2(Kab,Kcd,p,q,alpha,RPQ2,4+order,Boys_Table);
    g_SSSS_N = g_SSSS_N(order+1:4+order+1);
    end
    %The previous condition should increment the order of all functions,
    %without complicating the expressions, although _0 will no longer mean
    %_0
    
    gPSSS_0 = RPA*g_SSSS_N(1)+RWP*g_SSSS_N(2);
    gPSSS_1 = RPA*g_SSSS_N(2)+RWP*g_SSSS_N(3);
    gPSSS_2 = RPA*g_SSSS_N(3)+RWP*g_SSSS_N(4);
    gPSSS_3 = RPA*g_SSSS_N(4)+RWP*g_SSSS_N(5);
    
    gDm2SSS_0 = 1/2/p*(g_SSSS_N(1)-q/(p+q)*g_SSSS_N(2));
    gDm2SSS_1 = 1/2/p*(g_SSSS_N(2)-q/(p+q)*g_SSSS_N(3));
    gDm2SSS_2 = 1/2/p*(g_SSSS_N(3)-q/(p+q)*g_SSSS_N(4));
    
    
%     gDSSS_0 = [RPA(1)*gPSSS_0;RPA(2)*gPSSS_0(end-L1+2:end);RPA(3)*gPSSS_0(end)]...
%                 +[RWP(1)*gPSSS_1;RWP(2)*gPSSS_1(end-L1+2:end);RWP(3)*gPSSS_1(end)];
%      gDSSS_0 = [RPA*gPSSS_0(1);RPA(end-L1+2:end)*gPSSS_0(2);RPA(end)*gPSSS_0(end)]...
%                  +[RWP*gPSSS_1(1);RWP(end-L1+2:end)*gPSSS_1(2);RWP(end)*gPSSS_1(end)];
%             
%     gDSSS_0(1) = gDSSS_0(1) + gDm2SSS_0(1);       
%     gDSSS_0(4) = gDSSS_0(4) + gDm2SSS_0(2);
%     gDSSS_0(6) = gDSSS_0(6) + gDm2SSS_0(3);
%     
%     gDSSS_1 = [RPA(1)*gPSSS_1;RPA(2)*gPSSS_1(end-L1+2:end);RPA(3)*gPSSS_1(end)]...
%                 +[RWP(1)*gPSSS_2;RWP(2)*gPSSS_2(end-L1+2:end);RWP(3)*gPSSS_2(end)];
%      gDSSS_1 = [RPA*gPSSS_1(1);RPA(end-L1+2:end)*gPSSS_1(2);RPA(end)*gPSSS_1(end)]...
%                  +[RWP*gPSSS_2(1);RWP(end-L1+2:end)*gPSSS_2(2);RWP(end)*gPSSS_2(end)];
%             
%     gDSSS_1(1) = gDSSS_1(1) + gDm2SSS_1(1);       
%     gDSSS_1(4) = gDSSS_1(4) + gDm2SSS_1(2);
%     gDSSS_1(6) = gDSSS_1(6) + gDm2SSS_1(3);
%     
%      gDSSS_2 = [RPA(1)*gPSSS_2;RPA(2)*gPSSS_2(end-L1+2:end);RPA(3)*gPSSS_2(end)]...
%                  +[RWP(1)*gPSSS_3;RWP(2)*gPSSS_3(end-L1+2:end);RWP(3)*gPSSS_3(end)];
%      gDSSS_2 = [RPA*gPSSS_2(1);RPA(end-L1+2:end)*gPSSS_2(2);RPA(end)*gPSSS_2(end)]...
%                  +[RWP*gPSSS_3(1);RWP(end-L1+2:end)*gPSSS_3(2);RWP(end)*gPSSS_3(end)];
%             
%     gDSSS_2(1) = gDSSS_2(1) + gDm2SSS_2(1);       
%     gDSSS_2(4) = gDSSS_2(4) + gDm2SSS_2(2);
%     gDSSS_2(6) = gDSSS_2(6) + gDm2SSS_2(3);
    
     gDSSS_0 = [RPA(1)*gPSSS_0(1)+RWP(1)*gPSSS_1(1)+gDm2SSS_0
                RPA(2)*gPSSS_0(1)+RWP(2)*gPSSS_1(1)
                RPA(3)*gPSSS_0(1)+RWP(3)*gPSSS_1(1)
                RPA(2)*gPSSS_0(2)+RWP(2)*gPSSS_1(2)+gDm2SSS_0
                RPA(3)*gPSSS_0(2)+RWP(3)*gPSSS_1(2)
                RPA(3)*gPSSS_0(3)+RWP(3)*gPSSS_1(3)+gDm2SSS_0];

     gDSSS_1 = [RPA(1)*gPSSS_1(1)+RWP(1)*gPSSS_2(1)+gDm2SSS_1
                RPA(2)*gPSSS_1(1)+RWP(2)*gPSSS_2(1)
                RPA(3)*gPSSS_1(1)+RWP(3)*gPSSS_2(1)
                RPA(2)*gPSSS_1(2)+RWP(2)*gPSSS_2(2)+gDm2SSS_1
                RPA(3)*gPSSS_1(2)+RWP(3)*gPSSS_2(2)
                RPA(3)*gPSSS_1(3)+RWP(3)*gPSSS_2(3)+gDm2SSS_1];
            
     gDSSS_2 = [RPA(1)*gPSSS_2(1)+RWP(1)*gPSSS_3(1)+gDm2SSS_2
                RPA(2)*gPSSS_2(1)+RWP(2)*gPSSS_3(1)
                RPA(3)*gPSSS_2(1)+RWP(3)*gPSSS_3(1)
                RPA(2)*gPSSS_2(2)+RWP(2)*gPSSS_3(2)+gDm2SSS_2
                RPA(3)*gPSSS_2(2)+RWP(3)*gPSSS_3(2)
                RPA(3)*gPSSS_2(3)+RWP(3)*gPSSS_3(3)+gDm2SSS_2];        
    
    
%     ExpgLm2SSS_N = [gLm2SSS_N;gLm2SSS_N(end-(L1-1)+1:end);gLm2SSS_N(end)];
%     TwiceExpgLm2SSS_N = [ExpgLm2SSS_N;ExpgLm2SSS_N(end-L1+1:end);ExpgLm2SSS_N(end)];
    
    %According to Pritchard's SIMINT
    %3x1 vector
    gLm2SSS_N = 1/2/p*(gPSSS_0-q/(p+q)*gPSSS_1);
    gFSSS_0 =  [RPA(1)*gDSSS_0(1)+RWP(1)*gDSSS_1(1)+2*gLm2SSS_N(1)
                RPA(2)*gDSSS_0(1)+RWP(2)*gDSSS_1(1)
                RPA(3)*gDSSS_0(1)+RWP(3)*gDSSS_1(1)
                RPA(1)*gDSSS_0(4)+RWP(1)*gDSSS_1(4)
                RPA(3)*gDSSS_0(2)+RWP(3)*gDSSS_1(2)
                RPA(1)*gDSSS_0(6)+RWP(1)*gDSSS_1(6)
                RPA(2)*gDSSS_0(4)+RWP(2)*gDSSS_1(4)+2*gLm2SSS_N(2)
                RPA(3)*gDSSS_0(4)+RWP(3)*gDSSS_1(4)
                RPA(2)*gDSSS_0(6)+RWP(2)*gDSSS_1(6)
                RPA(3)*gDSSS_0(6)+RWP(3)*gDSSS_1(6)+2*gLm2SSS_N(3)];
    
    %3x1 vector
    gLm2SSS_Np1 = 1/2/p*(gPSSS_1-q/(p+q)*gPSSS_2);        
    gFSSS_1 =  [RPA(1)*gDSSS_1(1)+RWP(1)*gDSSS_2(1)+2*gLm2SSS_Np1(1)
                RPA(2)*gDSSS_1(1)+RWP(2)*gDSSS_2(1)
                RPA(3)*gDSSS_1(1)+RWP(3)*gDSSS_2(1)
                RPA(1)*gDSSS_1(4)+RWP(1)*gDSSS_2(4)
                RPA(3)*gDSSS_1(2)+RWP(3)*gDSSS_2(2)
                RPA(1)*gDSSS_1(6)+RWP(1)*gDSSS_2(6)
                RPA(2)*gDSSS_1(4)+RWP(2)*gDSSS_2(4)+2*gLm2SSS_Np1(2)
                RPA(3)*gDSSS_1(4)+RWP(3)*gDSSS_2(4)
                RPA(2)*gDSSS_1(6)+RWP(2)*gDSSS_2(6)
                RPA(3)*gDSSS_1(6)+RWP(3)*gDSSS_2(6)+2*gLm2SSS_Np1(3)];        
    
%     %I try a different ordering of the expansion.
%     TwiceExpgLm2SSS_N = [gLm2SSS_N(1)*2
%                         gLm2SSS_N(1)*0   
%                         gLm2SSS_N(1)*0
%                         gLm2SSS_N(1)*1
%                         gLm2SSS_N(1)*0
%                         gLm2SSS_N(1)*1
%                         gLm2SSS_N(2)*2
%                         gLm2SSS_N(2)*0
%                         gLm2SSS_N(2)*1
%                         gLm2SSS_N(3)*2];
%     gLSSS_N = [RPA(1)*gDSSS_0;RPA(2)*gDSSS_0(end-L1+1:end);RPA(3)*gDSSS_0(end)]...
%                 +[RWP(1)*gDSSS_1;RWP(2)*gDSSS_1(end-L1+1:end);RWP(3)*gDSSS_1(end)]...
%                 +TwiceExpgLm2SSS_N;%.*nz{L1};
%     term1N =   [RPA(1)*gDSSS_0(1)
%                 RPA(2)*gDSSS_0(1)
%                 RPA(3)*gDSSS_0(1)
%                 RPA(2)*gDSSS_0(2)
%                 RPA(3)*gDSSS_0(2)
%                 RPA(3)*gDSSS_0(3)
%                 RPA(2)*gDSSS_0(4)
%                 RPA(3)*gDSSS_0(4)
%                 RPA(3)*gDSSS_0(5)
%                 RPA(3)*gDSSS_0(6)];
%     term2N =   [RWP(1)*gDSSS_1(1)
%                 RWP(2)*gDSSS_1(1)
%                 RWP(3)*gDSSS_1(1)
%                 RWP(2)*gDSSS_1(2)
%                 RWP(3)*gDSSS_1(2)
%                 RWP(3)*gDSSS_1(3)
%                 RWP(2)*gDSSS_1(4)
%                 RWP(3)*gDSSS_1(4)
%                 RWP(3)*gDSSS_1(5)
%                 RWP(3)*gDSSS_1(6)];        
%     gLSSS_N = term1N+term2N+TwiceExpgLm2SSS_N;%.*nz{L1};
%             
%     gLm2SSS_Np1 = 1/2/p*(gPSSS_1-q/(p+q)*gPSSS_2);
%     ExpgLm2SSS_Np1 = [gLm2SSS_Np1;gLm2SSS_Np1(end-(L1-1)+1:end);gLm2SSS_Np1(end)];
%     TwiceExpgLm2SSS_Np1 = [ExpgLm2SSS_Np1;ExpgLm2SSS_Np1(end-L1+1:end);ExpgLm2SSS_Np1(end)];
% 
%         TwiceExpgLm2SSS_Np1 = [gLm2SSS_Np1(1)*2
%                                gLm2SSS_Np1(1)*0   
%                                gLm2SSS_Np1(1)*0
%                                gLm2SSS_Np1(1)*1
%                                gLm2SSS_Np1(1)*0
%                                gLm2SSS_Np1(1)*1
%                                gLm2SSS_Np1(2)*2
%                                gLm2SSS_Np1(2)*0
%                                gLm2SSS_Np1(2)*1
%                                gLm2SSS_Np1(3)*2];
%     
%     gLSSS_Np1 = [RPA(1)*gDSSS_1;RPA(2)*gDSSS_1(end-L1+1:end);RPA(3)*gDSSS_1(end)]...
%                 +[RWP(1)*gDSSS_2;RWP(2)*gDSSS_2(end-L1+1:end);RWP(3)*gDSSS_2(end)]...
%                 +TwiceExpgLm2SSS_Np1;%.*nz{L1};
%             
%     term1Np1 =   [RPA(1)*gDSSS_1(1)
%                 RPA(2)*gDSSS_1(1)
%                 RPA(3)*gDSSS_1(1)
%                 RPA(2)*gDSSS_1(2)
%                 RPA(3)*gDSSS_1(2)
%                 RPA(3)*gDSSS_1(3)
%                 RPA(2)*gDSSS_1(4)
%                 RPA(3)*gDSSS_1(4)
%                 RPA(3)*gDSSS_1(5)
%                 RPA(3)*gDSSS_1(6)];
%     term2Np1 =   [RWP(1)*gDSSS_2(1)
%                 RWP(2)*gDSSS_2(1)
%                 RWP(3)*gDSSS_2(1)
%                 RWP(2)*gDSSS_2(2)
%                 RWP(3)*gDSSS_2(2)
%                 RWP(3)*gDSSS_2(3)
%                 RWP(2)*gDSSS_2(4)
%                 RWP(3)*gDSSS_2(4)
%                 RWP(3)*gDSSS_2(5)
%                 RWP(3)*gDSSS_2(6)];        
%     gLSSS_Np1 = term1Np1+term2Np1+TwiceExpgLm2SSS_Np1;%.*nz{L1};        
    
    gLSSS_N = gFSSS_0;
    gLSSS_Np1 = gFSSS_1;
    gLm1SSS_N = gDSSS_0;
    gLm1SSS_Np1 = gDSSS_1;
            
elseif L1 == 4 %L1 > 3
    
    %[gLm1SSS_N,gLm1SSS_Np1,gLm2SSS_N,gLm2SSS_Np1] = vrr_1(a,b,c,d,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,3,0,0,0,Lmax,0,Boys_Table,nz);
    %[~,gLm1SSS_Np2,~,gLm2SSS_Np2] = vrr_1(a,b,c,d,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1-1,L2,L3,L4,Lmax,order+1,Boys_Table,nz);
    if order == 0
    g_SSSS_N = g_SSSS_N_2(Kab,Kcd,p,q,alpha,RPQ2,5,Boys_Table); %contains orders 0,1,2,3
    else
    g_SSSS_N = g_SSSS_N_2(Kab,Kcd,p,q,alpha,RPQ2,5+order,Boys_Table);
    g_SSSS_N = g_SSSS_N(order+1:5+order+1);
    end
    
    gPSSS_0 = RPA*g_SSSS_N(1)+RWP*g_SSSS_N(2);
    gPSSS_1 = RPA*g_SSSS_N(2)+RWP*g_SSSS_N(3);
    gPSSS_2 = RPA*g_SSSS_N(3)+RWP*g_SSSS_N(4);
    gPSSS_3 = RPA*g_SSSS_N(4)+RWP*g_SSSS_N(5);
    gPSSS_4 = RPA*g_SSSS_N(5)+RWP*g_SSSS_N(6);
    
    gDm2SSS_0 = 1/2/p*(g_SSSS_N(1)-q/(p+q)*g_SSSS_N(2));
    gDm2SSS_1 = 1/2/p*(g_SSSS_N(2)-q/(p+q)*g_SSSS_N(3));
    gDm2SSS_2 = 1/2/p*(g_SSSS_N(3)-q/(p+q)*g_SSSS_N(4));
    gDm2SSS_3 = 1/2/p*(g_SSSS_N(4)-q/(p+q)*g_SSSS_N(5));
    
     gDSSS_0 = [RPA(1)*gPSSS_0(1)+RWP(1)*gPSSS_1(1)+gDm2SSS_0
                RPA(2)*gPSSS_0(1)+RWP(2)*gPSSS_1(1)
                RPA(3)*gPSSS_0(1)+RWP(3)*gPSSS_1(1)
                RPA(2)*gPSSS_0(2)+RWP(2)*gPSSS_1(2)+gDm2SSS_0
                RPA(3)*gPSSS_0(2)+RWP(3)*gPSSS_1(2)
                RPA(3)*gPSSS_0(3)+RWP(3)*gPSSS_1(3)+gDm2SSS_0];

     gDSSS_1 = [RPA(1)*gPSSS_1(1)+RWP(1)*gPSSS_2(1)+gDm2SSS_1
                RPA(2)*gPSSS_1(1)+RWP(2)*gPSSS_2(1)
                RPA(3)*gPSSS_1(1)+RWP(3)*gPSSS_2(1)
                RPA(2)*gPSSS_1(2)+RWP(2)*gPSSS_2(2)+gDm2SSS_1
                RPA(3)*gPSSS_1(2)+RWP(3)*gPSSS_2(2)
                RPA(3)*gPSSS_1(3)+RWP(3)*gPSSS_2(3)+gDm2SSS_1];
            
     gDSSS_2 = [RPA(1)*gPSSS_2(1)+RWP(1)*gPSSS_3(1)+gDm2SSS_2
                RPA(2)*gPSSS_2(1)+RWP(2)*gPSSS_3(1)
                RPA(3)*gPSSS_2(1)+RWP(3)*gPSSS_3(1)
                RPA(2)*gPSSS_2(2)+RWP(2)*gPSSS_3(2)+gDm2SSS_2
                RPA(3)*gPSSS_2(2)+RWP(3)*gPSSS_3(2)
                RPA(3)*gPSSS_2(3)+RWP(3)*gPSSS_3(3)+gDm2SSS_2];
            
     gDSSS_3 = [RPA(1)*gPSSS_3(1)+RWP(1)*gPSSS_4(1)+gDm2SSS_3
                RPA(2)*gPSSS_3(1)+RWP(2)*gPSSS_4(1)
                RPA(3)*gPSSS_3(1)+RWP(3)*gPSSS_4(1)
                RPA(2)*gPSSS_3(2)+RWP(2)*gPSSS_4(2)+gDm2SSS_3
                RPA(3)*gPSSS_3(2)+RWP(3)*gPSSS_4(2)
                RPA(3)*gPSSS_3(3)+RWP(3)*gPSSS_4(3)+gDm2SSS_3];      

%{
%     DiffgLm2SSS_N = 1/2/p*(gLm2SSS_N-q/(p+q)*gLm2SSS_Np1);
%     ExpgLm2SSS_N = [DiffgLm2SSS_N;DiffgLm2SSS_N(end-(L1-1)+1:end);DiffgLm2SSS_N(end)];
%     TwiceExpgLm2SSS_N = [ExpgLm2SSS_N;ExpgLm2SSS_N(end-L1+1:end);ExpgLm2SSS_N(end)];
%     
%     TwiceExpgLm2SSS_N = [DiffgLm2SSS_N(1)*3
%                         DiffgLm2SSS_N(1)*0
%                         DiffgLm2SSS_N(1)*0
%                         DiffgLm2SSS_N(1)*1
%                         DiffgLm2SSS_N(1)*0
%                         DiffgLm2SSS_N(1)*1
%                         DiffgLm2SSS_N(2)*2
%                         DiffgLm2SSS_N(2)*0
%                         DiffgLm2SSS_N(2)*1
%                         DiffgLm2SSS_N(3)*2
%                         DiffgLm2SSS_N(4)*3
%                         DiffgLm2SSS_N(4)*0
%                         DiffgLm2SSS_N(4)*1
%                         DiffgLm2SSS_N(5)*2
%                         DiffgLm2SSS_N(6)*3];
%     
%     DiffgLm2SSS_Np1 = 1/2/p*(gLm2SSS_Np1-q/(p+q)*gLm2SSS_Np2);
%     ExpgLm2SSS_Np1 = [DiffgLm2SSS_Np1;DiffgLm2SSS_Np1(end-(L1-1)+1:end);DiffgLm2SSS_Np1(end)];
%     TwiceExpgLm2SSS_Np1 = [ExpgLm2SSS_Np1;ExpgLm2SSS_Np1(end-L1+1:end);ExpgLm2SSS_Np1(end)];
%      
%     TwiceExpgLm2SSS_Np1 = [DiffgLm2SSS_Np1(1)*3
%                         DiffgLm2SSS_Np1(1)*0
%                         DiffgLm2SSS_Np1(1)*0
%                         DiffgLm2SSS_Np1(1)*1
%                         DiffgLm2SSS_Np1(1)*0
%                         DiffgLm2SSS_Np1(1)*1
%                         DiffgLm2SSS_Np1(2)*2
%                         DiffgLm2SSS_Np1(2)*0
%                         DiffgLm2SSS_Np1(2)*1
%                         DiffgLm2SSS_Np1(3)*2
%                         DiffgLm2SSS_Np1(4)*3
%                         DiffgLm2SSS_Np1(4)*0
%                         DiffgLm2SSS_Np1(4)*1
%                         DiffgLm2SSS_Np1(5)*2
%                         DiffgLm2SSS_Np1(6)*3];
%     
% %     gLSSS_N = [RPA(1)*gLm1SSS_N;RPA(2)*gLm1SSS_N(end-L1+1:end);RPA(3)*gLm1SSS_N(end)]...
% %                 +[RWP(1)*gLm1SSS_Np1;RWP(2)*gLm1SSS_Np1(end-L1+1:end);RWP(3)*gLm1SSS_Np1(end)]...
% %                 +TwiceExpgLm2SSS_N;%.*nz{L1};
%             
%     term1N =   [RPA(1)*gLm1SSS_N(1)
%                 RPA(2)*gLm1SSS_N(1)
%                 RPA(3)*gLm1SSS_N(1)
%                 RPA(2)*gLm1SSS_N(2)
%                 RPA(3)*gLm1SSS_N(2)
%                 RPA(3)*gLm1SSS_N(3)
%                 RPA(2)*gLm1SSS_N(4)
%                 RPA(3)*gLm1SSS_N(4)
%                 RPA(3)*gLm1SSS_N(5)
%                 RPA(3)*gLm1SSS_N(6)
%                 RPA(2)*gLm1SSS_N(7)
%                 RPA(3)*gLm1SSS_N(7)
%                 RPA(3)*gLm1SSS_N(8)
%                 RPA(3)*gLm1SSS_N(9)
%                 RPA(3)*gLm1SSS_N(10)];
%     term2N =   [RWP(1)*gLm1SSS_Np1(1)
%                 RWP(2)*gLm1SSS_Np1(1)
%                 RWP(3)*gLm1SSS_Np1(1)
%                 RWP(2)*gLm1SSS_Np1(2)
%                 RWP(3)*gLm1SSS_Np1(2)
%                 RWP(3)*gLm1SSS_Np1(3)
%                 RWP(2)*gLm1SSS_Np1(4)
%                 RWP(3)*gLm1SSS_Np1(4)
%                 RWP(3)*gLm1SSS_Np1(5)
%                 RWP(3)*gLm1SSS_Np1(6)
%                 RWP(2)*gLm1SSS_Np1(7)
%                 RWP(3)*gLm1SSS_Np1(7)
%                 RWP(3)*gLm1SSS_Np1(8)
%                 RWP(3)*gLm1SSS_Np1(9)
%                 RWP(3)*gLm1SSS_Np1(10)];        
% %     gLSSS_N = [RPA(1)*gLm1SSS_N(1);RPA(end-L1+2:end)*gLm1SSS_N(2);RPA(end)*gLm1SSS_N(end)]...
% %                  +[RWP*gLm1SSS_Np1(1);RWP(end-L1+2:end)*gLm1SSS_Np1(2);RWP(end)*gLm1SSS_Np1(end)]...
% %                  +TwiceExpgLm2SSS_N;%.*nz{L1};          
%     gLSSS_N = term1N+term2N+TwiceExpgLm2SSS_N;%.*nz{L1};
% %     gLSSS_Np1 = [RPA(1)*gLm1SSS_Np1(1);RPA(end-L1+2:end)*gLm1SSS_Np1(2);RPA(end)*gLm1SSS_Np1(end)]...
% %                  +[RWP*gLm1SSS_Np2(1);RWP(end-L1+2:end)*gLm1SSS_Np2(2);RWP(end)*gLm1SSS_Np2(end)]...
% %                  +TwiceExpgLm2SSS_Np1;%.*nz{L1};
%}
    gFm2SSS_N = 1/2/p*(gPSSS_0-q/(p+q)*gPSSS_1);
    gFSSS_0 =  [RPA(1)*gDSSS_0(1)+RWP(1)*gDSSS_1(1)+2*gFm2SSS_N(1)
                RPA(2)*gDSSS_0(1)+RWP(2)*gDSSS_1(1)
                RPA(3)*gDSSS_0(1)+RWP(3)*gDSSS_1(1)
                RPA(1)*gDSSS_0(4)+RWP(1)*gDSSS_1(4)
                RPA(3)*gDSSS_0(2)+RWP(3)*gDSSS_1(2)
                RPA(1)*gDSSS_0(6)+RWP(1)*gDSSS_1(6)
                RPA(2)*gDSSS_0(4)+RWP(2)*gDSSS_1(4)+2*gFm2SSS_N(2)
                RPA(3)*gDSSS_0(4)+RWP(3)*gDSSS_1(4)
                RPA(2)*gDSSS_0(6)+RWP(2)*gDSSS_1(6)
                RPA(3)*gDSSS_0(6)+RWP(3)*gDSSS_1(6)+2*gFm2SSS_N(3)];
    
    %3x1 vector
    gFm2SSS_Np1 = 1/2/p*(gPSSS_1-q/(p+q)*gPSSS_2);        
    gFSSS_1 =  [RPA(1)*gDSSS_1(1)+RWP(1)*gDSSS_2(1)+2*gFm2SSS_Np1(1)
                RPA(2)*gDSSS_1(1)+RWP(2)*gDSSS_2(1)
                RPA(3)*gDSSS_1(1)+RWP(3)*gDSSS_2(1)
                RPA(1)*gDSSS_1(4)+RWP(1)*gDSSS_2(4)
                RPA(3)*gDSSS_1(2)+RWP(3)*gDSSS_2(2)
                RPA(1)*gDSSS_1(6)+RWP(1)*gDSSS_2(6)
                RPA(2)*gDSSS_1(4)+RWP(2)*gDSSS_2(4)+2*gFm2SSS_Np1(2)
                RPA(3)*gDSSS_1(4)+RWP(3)*gDSSS_2(4)
                RPA(2)*gDSSS_1(6)+RWP(2)*gDSSS_2(6)
                RPA(3)*gDSSS_1(6)+RWP(3)*gDSSS_2(6)+2*gFm2SSS_Np1(3)];   
      
    gFm2SSS_Np2 = 1/2/p*(gPSSS_2-q/(p+q)*gPSSS_3);        
    gFSSS_2 =  [RPA(1)*gDSSS_2(1)+RWP(1)*gDSSS_3(1)+2*gFm2SSS_Np2(1)
                RPA(2)*gDSSS_2(1)+RWP(2)*gDSSS_3(1)
                RPA(3)*gDSSS_2(1)+RWP(3)*gDSSS_3(1)
                RPA(1)*gDSSS_2(4)+RWP(1)*gDSSS_3(4)
                RPA(3)*gDSSS_2(2)+RWP(3)*gDSSS_3(2)
                RPA(1)*gDSSS_2(6)+RWP(1)*gDSSS_3(6)
                RPA(2)*gDSSS_2(4)+RWP(2)*gDSSS_3(4)+2*gFm2SSS_Np2(2)
                RPA(3)*gDSSS_2(4)+RWP(3)*gDSSS_3(4)
                RPA(2)*gDSSS_2(6)+RWP(2)*gDSSS_3(6)
                RPA(3)*gDSSS_2(6)+RWP(3)*gDSSS_3(6)+2*gFm2SSS_Np2(3)];        
                
                %6 x 1 vector     
     
     
     gLm1SSS_N = gFSSS_0;
     gLm1SSS_Np1 = gFSSS_1;
     
     DiffgLm2SSS_N = 1/2/p*(gDSSS_0-q/(p+q)*gDSSS_1);
     gLSSS_N = [RPA(1)*gLm1SSS_N(1)+RWP(1)*gLm1SSS_Np1(1)+3*DiffgLm2SSS_N(1)
                RPA(2)*gLm1SSS_N(1)+RWP(2)*gLm1SSS_Np1(1)
                RPA(3)*gLm1SSS_N(1)+RWP(3)*gLm1SSS_Np1(1)
                RPA(2)*gLm1SSS_N(2)+RWP(2)*gLm1SSS_Np1(2)+1*DiffgLm2SSS_N(1)
                RPA(3)*gLm1SSS_N(2)+RWP(3)*gLm1SSS_Np1(2)
                RPA(3)*gLm1SSS_N(3)+RWP(3)*gLm1SSS_Np1(3)+1*DiffgLm2SSS_N(1)
                RPA(1)*gLm1SSS_N(7)+RWP(1)*gLm1SSS_Np1(7)
                RPA(3)*gLm1SSS_N(4)+RWP(3)*gLm1SSS_Np1(4)
                RPA(2)*gLm1SSS_N(6)+RWP(2)*gLm1SSS_Np1(6)
                RPA(1)*gLm1SSS_N(10)+RWP(1)*gLm1SSS_Np1(10)
                RPA(2)*gLm1SSS_N(7)+RWP(2)*gLm1SSS_Np1(7)+3*DiffgLm2SSS_N(4)
                RPA(3)*gLm1SSS_N(7)+RWP(3)*gLm1SSS_Np1(7)
                RPA(3)*gLm1SSS_N(8)+RWP(3)*gLm1SSS_Np1(8)+1*DiffgLm2SSS_N(4)
                RPA(2)*gLm1SSS_N(10)+RWP(2)*gLm1SSS_Np1(10)
                RPA(3)*gLm1SSS_N(10)+RWP(3)*gLm1SSS_Np1(10)+3*DiffgLm2SSS_N(6)];

    gLSSS_Np1 = [];
else %L1 > 4
    gLSSS_N = zeros((L1+1)*(L1+2)/2,1);
    gLSSS_Np1 = zeros((L1+1)*(L1+2)/2,1);
    gLm1SSS_N = zeros((L1)*(L1+1)/2,1);
    gLm1SSS_Np1 = zeros((L1)*(L1+1)/2,1);
end

end
