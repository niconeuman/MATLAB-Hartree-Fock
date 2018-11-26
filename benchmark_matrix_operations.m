close all


Size_1 = 100;
Size_2 = 400;
Size_3 = 2000;
Sizes = [(1:Size_1) (Size_1+1:3:Size_2) (Size_2+1:5:Size_3)]';%'
Nelements =  length(Sizes);
time_rand_vector_gen = zeros(Nelements,1);
time_outer_product = zeros(Nelements,1);
time_Sum = zeros(Nelements,1);
time_mult_add = zeros(Nelements,1);
time_ps_int = zeros(Nelements,1);
time_pp_int = zeros(Nelements,1);




for el = 1:Nelements
        Size = Sizes(el);
        if Size < 1000
                Nrep = 100;
        elseif Size < 1500
                Nrep = 20;
        else
                Nrep = 5;
        end
                tcnt = 0;
                for k = 1:Nrep
                        tic;
                        vector1 = rand(Size,1);
                        vector2 = rand(Size,1);
                        tcnt = tcnt + toc;
                end

                time_rand_vector_gen(el) = tcnt/Nrep;

                tcnt = 0;
                        for k = 1:Nrep
                                vector1 = rand(Size,1);
                                vector2 = rand(Size,1);

                                tic;
                                Outer_product = vector1*vector2'; %'
                                tcnt = tcnt + toc;
                        end

                time_outer_product(el) =  tcnt/Nrep;

                tcnt = 0;
                        for k = 1:Nrep
                                vector1 = rand(Size,1);
                                vector2 = rand(Size,1);
                                tic;
                                Sum = vector1+vector2;
                                tcnt = tcnt + toc;
                        end

                time_Sum(el) =  tcnt/Nrep;

                tcnt = 0;
                        for k = 1:Nrep
                                vector1 = rand(Size,1);
                                vector2 = rand(Size,1);
                                tic;
                                Result = vector1+vector1.*vector2;
                                tcnt = tcnt + toc;
                        end

                time_mult_add(el) =  tcnt/Nrep;

                tcnt = 0;
                        for k = 1:Nrep
                                RPA = rand(Size,1);
                                gSSSS0 = rand(Size,1);
                                RWP = rand(Size,1);
                                gSSSS1 = rand(Size,1);
                                tic;
                                Result = RPA.*gSSSS0+RWP.*gSSSS1;
                                Result2 = sum(Result);
                                tcnt = tcnt + 3*toc; %3 is because there are 3 components
                        end

                time_ps_int(el) =  tcnt/Nrep;

                tcnt = 0;
                        for k = 1:Nrep
                                RPAx = rand(Size,1);
                                RPBx = rand(Size,1);
                                RPAy = rand(Size,1);
                                RPBy = rand(Size,1);
                                RPAz = rand(Size,1);
                                RPBz = rand(Size,1);
                                gSSSS0 = rand(Size,1);
                                RWPx = rand(Size,1);
                                RWPy = rand(Size,1);
                                RWPz = rand(Size,1);
                                gSSSS1 = rand(Size,1);
                                pValues = rand(Size,1);
                                qValues = rand(Size,1);
                                ppqValues = rand(Size,1);
                                gSSSS2 = rand(Size,1);
                                tic;
                                pxsss0 = RPAx.*gSSSS0+RWPx.*gSSSS1;
                                pxsss1 = RPAx.*gSSSS1+RWPx.*gSSSS2;
                                pxpxss0 = RPBx.*pxsss0+RWPx.*pxsss1+0.5./pValues.**(gSSSS0-qValues./ppqValues.*gSSSS1);
                                pxpyss0 = RPBy.*pxsss0+RWPy.*pxsss1;
                                pxpzss0 = RPBz.*pxsss0+RWPz.*pxsss1;

                                Resultx = sum(pxpxss0);
                                Resulty = sum(pxpyss0);
                                Resultz = sum(pxpzss0);
                                tcnt = tcnt + 3*toc; %3 is because there are 3 components
                        end

                time_pp_int(el) =  tcnt/Nrep;

end
mod_sizes_2 = mod(Sizes,2);
mod_sizes_4 = mod(Sizes,4);


plot(Sizes,time_outer_product,'-dk','MarkerFaceColor','k');
plot(Sizes(mod_sizes_4==0),time_outer_product(mod_sizes_4==0),'-dr','MarkerFaceColor','r'); hold on
plot(Sizes(mod_sizes_2==1),time_outer_product(mod_sizes_2==1),'-db','MarkerFaceColor','b');
legend('Vector Outer Product');
xlabel('Vector Size');
ylabel('Time (s)');
figure
plot(Sizes,time_Sum,'-dk','MarkerFaceColor','k');
plot(Sizes(mod_sizes_4==0),time_Sum(mod_sizes_4==0),'-dr','MarkerFaceColor','r'); hold on
plot(Sizes(mod_sizes_2==1),time_Sum(mod_sizes_2==1),'-db','MarkerFaceColor','b');
legend('Vector Sum');
xlabel('Vector Size');
ylabel('Time (s)');
figure
plot(Sizes,time_mult_add,'-dk','MarkerFaceColor','k'); hold on
plot(Sizes(mod_sizes_4==0),time_mult_add(mod_sizes_4==0),'-dr','MarkerFaceColor','r'); hold on
plot(Sizes(mod_sizes_2==1),time_mult_add(mod_sizes_2==1),'-db','MarkerFaceColor','b');
legend('Vector Elementwise Product');
xlabel('Vector Size');
ylabel('Time (s)');
figure
plot(Sizes,time_ps_int,'-dk','MarkerFaceColor','k'); hold on
plot(Sizes(mod_sizes_4==0),time_ps_int(mod_sizes_4==0),'-dr','MarkerFaceColor','r'); hold on
plot(Sizes(mod_sizes_2==1),time_ps_int(mod_sizes_2==1),'-db','MarkerFaceColor','b');
legend('psss type integral');
xlabel('Vector Size');
ylabel('Time (s)');
figure
plot(Sizes,time_pp_int,'-dk','MarkerFaceColor','k'); hold on
plot(Sizes(mod_sizes_4==0),time_pp_int(mod_sizes_4==0),'-dr','MarkerFaceColor','r'); hold on
plot(Sizes(mod_sizes_2==1),time_pp_int(mod_sizes_2==1),'-db','MarkerFaceColor','b');
legend('ppss type integral');
xlabel('Vector Size');
ylabel('Time (s)');
