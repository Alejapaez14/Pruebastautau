library(shiny)
library(ggplot2)
library(tseries)
library(bslib)
library(htmltools)
library(thematic)
library(rsconnect)



ui<-fluidPage(
  deployApp(),
        theme = bs_theme_update(theme = bs_theme(fg = "rgb(11, 233, 122)", 
                                              primary = "#93BE06", 
                                              secondary = "#FADBD8", 
                                              base_font = font_collection("Rubik Bubbles","cursive"),
                                              code_font = font_collection(font_google("Fredericka the Great"),"cursive"), 
                                              heading_font = font_collection(font_google("Fredericka the Great"),"cursive"), 
                                              font_scale = NULL, bootswatch = "cyborg", 
                                              bg = "rgb(6, 6, 6)")),
              navbarPage(titlePanel("Prueba Tau~Tau, Tau~mu, Tau")),
              HTML("<body>Para esta Shinny App tendremos las tres puebas cada una estar?? en cada tab.
                    Se deber?? elegir el n??mero de escenarios y el n??mero de observaciones con los cuales se desea evaluar la prueba</body>"),
              # Sidebar with a slider input for number of bins
              tabsetPanel(tabPanel("Prueba Tau",
                                   sidebarLayout(
                                     sidebarPanel(
                                       tags$style(
                                         ".back{
                                         background: -ms-linear-gradient(left, rgba(255,255,255,0) 0%,rgba(255,255,255,0.8) 50%,rgba(128,186,232,0) 99%,rgba(125,185,232,0) 100%);
                                         }
                                         .green{
                                              color: #93BE06;
                                           }
                                            .red{
                                              color: #F72828;
                                            }
                                          .purple{
                                              color: #780EF5;
                                           }
                                           
                                            "),
                                       h3(p(class = "green", "Par??metros a elegir para la Prueba Tau ~ Mu")),
                                       numericInput("Al", p(class = "green", "N??mero m??ximo que desea que el slider tenga"), value = 20 ),
                                       sliderInput("random",p(class = "green", "N??mero de observaciones"), min = 0, max = 100000, value =2),
                                       sliderInput("esce",p(class = "green", "N??mero de escenarios"), min = 0, max = 100000, value =2)
                                     ),
                                     mainPanel(
                                       h3(p(class = "green","Resultados gr??ficos de la prueba Tau")),
                                       plotOutput("Tau"),
                                     )
                                   )
                                   
              ),
              tabPanel("Prueba Tau~Mu",
                       sidebarLayout(
                         sidebarPanel(
                           h3(p(class = "purple", "Par??metros a elegir para la Prueba Tau ~ Mu")),
                           numericInput("bl", p(class = "purple", "N??mero m??ximo que desea que el slider tenga"), value = 20),
                           sliderInput("Esce",p(class = "purple", "N??mero de escenarios"), min = 0, max = 100000, value =2),
                           sliderInput("Random",p(class = "purple", "N??mero de observaciones"), min = 0, max = 100000, value =2),
                         ),
                         mainPanel(
                           h3(p(class = "purple","Resultados gr??ficos de la prueba Tau~Mu")),
                           plotOutput("Taumu"),
                         )
                       )
              ),
              tabPanel("Prueba Tau~Tau",
                       sidebarLayout(
                         sidebarPanel(
                           h3(p(class = "red", "P??rametros a elegir para la Prueba Tau ~ Mu")),
                           numericInput("cl", p(class = "red", "N??mero m??ximo que desea que el slider tenga"), value = 20),
                           sliderInput("RANDOM",p(class = "red", "N??mero de observaciones"), min = 0, max = 100000, value =2),
                           sliderInput("ESCE",p(class = "red", "N??mero de escenarios"), min = 0, max = 100000, value =2),
                         ),
                         mainPanel(
                           h3(p(class = "red","Resultados gr??ficos de la prueba Tau~Tau")),
                           plotOutput("TauTau"),
                         )
                       )
              ),
              )             
)



server <- function(input,output,session){
  observeEvent(input$Al,{
    updateSliderInput(session, "random", max = input$Al)
    updateSliderInput(session, "esce", max = input$Al)
  })
  observeEvent(input$bl,{
    updateSliderInput(session, "Random", max = input$bl)
    updateSliderInput(session, "Esce", max = input$bl)
  })
  observeEvent(input$cl,{
    updateSliderInput(session, "RANDOM", max = input$cl)
    updateSliderInput(session, "ESCE", max = input$cl)
  })
  
  Tau = reactive({
    #Inputs
    T = input$random
    S = input$esce 
    Z_0 = 0
    Sigma_a = 1
    Zt <- c(Z_0)                                                    #Creando matriz con el valor inicial
    
    #Creando una matriz vac??a para almacenar los estad??sticos
    Tau = c()
    
    for (n in 1 : S){
      #Generando los errores
      set.seed(n)
      At <- rnorm(n = T ,mean = 0, sd = Sigma_a)
      
      #Corriendo la funci??n generadora del PGD
      for(t in 2 : T){
        Zt[t] = Zt[t - 1] + At[t]
      }
      
      #Calculando la serie diferenciada
      Diferencia = function(Data){                                 # Operador de diferencia
        Serie.dif = matrix(data = NA, nrow = (T - 1), ncol = 1)
        for(i in 1 : (T - 1)){
          Serie.dif[i] <- Data[i + 1] - Data[i]
        }
        return(Serie.dif)
      }
      X.Diff = Diferencia(Data = Zt)
      
      #Calculando la serie rezagada
      X.Rez <- Zt[1 : (T - 1)]
      
      #Estimando la regresi??n (OLS)
      Y = X.Diff
      X = cbind(X.Rez)
      Coeff = solve(t(X) %*% X) %*% (t(X) %*% Y)
      Gamma = Coeff[1]                                               # Coeficiente Gamma estimado
      
      # Estimaci??n sigma^2 (Beta mo??o de OLS)
      Ym = X %*% Coeff
      Em = Y - Ym
      Sigma2 = as.numeric((t(Em) %*% Em) / (T - 1))
      Var.cov = (Sigma2 * solve(t(X) %*% X))
      Gamma_std = sqrt(Var.cov[1, 1])
      
      #Llenando la matriz
      Tau[n] <- Gamma / Gamma_std
    }
    # Quantiles
    Q = matrix(data = c(quantile(Tau, probs = 0.2), quantile(Tau, probs = 0.1), quantile(Tau, probs = 0.05), quantile(Tau, probs = 0.01)), ncol = 4, nrow = 1)
    colnames(Q) = c('20%', '10%', '5%', '1%')
    
    # Gr??ficas
    # print(G_Tau)
    Den = density(Tau)
    value1 = Q[, 1]
    value2 = Q[, 2]
    value3 = Q[, 3]
    value4 = Q[, 4]
    plot(Den, main = sprintf('Dist. Tau con %s datos y %s escenarios', T, S), xlab = "Estad??stico Tau", ylab = "Densidad", col = "black", lty = 1, lwd = 4)
    legend(x= "topright" , legend = c("Alfa = 0.2", Q[, 1], "Alfa = 0.1", Q[, 2], "Alfa = 0.05", Q[, 3], "Alfa = 0.01", Q[, 4]))
    polygon(c(Den$x[Den$x <= value1], value1),
            c(Den$y[Den$x <= value1], 0),
            col = rgb(0, 1, 0, alpha = 0.2),
            border = 1)
    polygon(c(Den$x[Den$x <= value2], value2),
            c(Den$y[Den$x <= value2], 0),
            col = rgb(0, 1, 0, alpha = 0.2),
            border = 1)
    polygon(c(Den$x[Den$x <= value3], value3),
            c(Den$y[Den$x <= value3], 0),
            col = rgb(0, 1, 0, alpha = 0.2),
            border = 1)
    polygon(c(Den$x[Den$x <= value4], value4),
            c(Den$y[Den$x <= value4], 0),
            col = rgb(0, 1, 0, alpha = 0.2),
            border = 1)
    abline(v = mean(Tau), col = "green", lwd = 3, lty = 2)
    
  })
  Taumu =  reactive({
    #Input 
    T = input$Random
    S = input$Esce
    Z_0 = 0
    Sigma_a = 1
    Zt <- c(Z_0)                                                    #Creando matriz con el valor inicial
    
    #Creando una matriz vac??a para almacenar los estad??sticos
    Tau_M = c()
    
    for (n in 1 : S){
      #Generando los errores
      set.seed(n)
      At <- rnorm(n = T ,mean = 0, sd = Sigma_a)
      
      #Corriendo la funci??n generadora del PGD
      for(t in 2 : T){
        Zt[t] = Zt[t - 1] + At[t]
      }
      
      #Calculando la serie diferenciada
      Diferencia = function(Data){                                 # Operador de diferencia
        Serie.dif = matrix(data = NA, nrow = (T - 1), ncol = 1)
        for(i in 1 : (T - 1)){
          Serie.dif[i] <- Data[i + 1] - Data[i]
        }
        return(Serie.dif)
      }
      X.Diff = Diferencia(Data = Zt)
      
      #Calculando la serie rezagada
      X.Rez <- Zt[1 : (T - 1)]
      
      #Estimando la regresi??n (OLS)
      Y = X.Diff
      X = cbind(1, X.Rez)
      Coeff = solve(t(X) %*% X) %*% (t(X) %*% Y)
      Gamma = Coeff[2]                                               # Coeficiente Gamma estimado
      
      # Estimaci??n sigma^2 (Beta mo??o de OLS)
      Ym = X %*% Coeff
      Em = Y - Ym
      Sigma2 = as.numeric((t(Em) %*% Em) / (T - 2))
      Var.cov = (Sigma2 * solve(t(X) %*% X))
      Gamma_std = sqrt(Var.cov[2, 2])
      
      #Llenando la matriz
      Tau_M[n] <- Gamma / Gamma_std
    }
    # Quantiles
    Q = matrix(data = c(quantile(Tau_M, probs = 0.2), quantile(Tau_M, probs = 0.1), quantile(Tau_M, probs = 0.05), quantile(Tau_M, probs = 0.01)), ncol = 4, nrow = 1)
    colnames(Q) = c('20%', '10%', '5%', '1%')
    
    # Gr??ficas
    # print(G_Tau_M)
    Den = density(Tau_M)
    value1 = Q[, 1]
    value2 = Q[, 2]
    value3 = Q[, 3]
    value4 = Q[, 4]
    plot(Den, main = sprintf('Dist. Tau ~ Mu con %s datos y %s escenarios', T, S), xlab="Estad??stico Tau - Mu", ylab="Densidad", col="black", lty=1, lwd=4)
    legend(x= "topright" , legend = c("Alfa = 0.2", Q[, 1], "Alfa = 0.1", Q[, 2], "Alfa = 0.05", Q[, 3], "Alfa = 0.01", Q[, 4]))
    polygon(c(Den$x[Den$x <= value1], value1),
            c(Den$y[Den$x <= value1], 0),
            col = rgb(0, 0, 1, alpha = 0.2),
            border = 1)
    polygon(c(Den$x[Den$x <= value2], value2),
            c(Den$y[Den$x <= value2], 0),
            col = rgb(0, 0, 1, alpha = 0.2),
            border = 1)
    polygon(c(Den$x[Den$x <= value3], value3),
            c(Den$y[Den$x <= value3], 0),
            col = rgb(0, 0, 1, alpha = 0.2),
            border = 1)
    polygon(c(Den$x[Den$x <= value4], value4),
            c(Den$y[Den$x <= value4], 0),
            col = rgb(0, 0, 1, alpha = 0.2),
            border = 1)
    abline(v = mean(Tau_M), col = "blue", lwd = 3, lty = 2)
  })
  TauTau = reactive({
    #Input 
    T = input$RANDOM
    S = input$ESCE
    Z_0 = 0
    Sigma_a = 1
    Zt <- c(Z_0)                                                    #Creando matriz con el valor inicial
    
    #Creando una matriz vac??a para almacenar los estad??sticos
    Tau_T = c()
    
    for (n in 1 : S){
      #Generando los errores
      set.seed(n)
      At <- rnorm(n = T ,mean = 0, sd = Sigma_a)
      
      #Corriendo la funci??n generadora del PGD
      for(t in 2 : T){
        Zt[t] = Zt[t - 1] + At[t]
      }
      
      #Calculando la serie diferenciada
      Diferencia = function(Data){                                 # Operador de diferencia
        Serie.dif = matrix(data = NA, nrow = (T - 1), ncol = 1)
        for(i in 1 : (T - 1)){
          Serie.dif[i] <- Data[i + 1] - Data[i]
        }
        return(Serie.dif)
      }
      X.Diff = Diferencia(Data = Zt)
      
      #Calculando la serie rezagada
      X.Rez <- Zt[1 : (T - 1)]
      
      #Estimando la regresi??n (OLS)
      Y = X.Diff
      X = cbind(1, c(1 : (T - 1)), X.Rez)
      Coeff = solve(t(X) %*% X) %*% (t(X) %*% Y)
      Gamma = Coeff[3]                                               # Coeficiente Gamma estimado
      
      # Estimaci??n sigma^2 (Beta mo??o de OLS)
      Ym = X %*% Coeff
      Em = Y - Ym
      Sigma2 = as.numeric((t(Em) %*% Em) / (T - 3))
      Var.cov = (Sigma2 * solve(t(X) %*% X))
      Gamma_std = sqrt(Var.cov[3, 3])
      
      #Llenando la matriz
      Tau_T[n] <- Gamma / Gamma_std
    }
    
    # Quantiles
    Q = matrix(data = c(quantile(Tau_T, probs = 0.2), quantile(Tau_T, probs = 0.1), quantile(Tau_T, probs = 0.05), quantile(Tau_T, probs = 0.01)), ncol = 4, nrow = 1)
    colnames(Q) = c('20%', '10%', '5%', '1%')
    
    # Gr??ficas
    # print(G_Tau_T)
    Den = density(Tau_T)
    value1 = Q[, 1]
    value2 = Q[, 2]
    value3 = Q[, 3]
    value4 = Q[, 4]
    plot(Den, main = sprintf('Dist. Tau ~ T con %s datos y %s escenarios', T, S), xlab = "Estad??stico Tau - Tau", ylab = "Densidad", col = "black", lty = 1, lwd = 4)
    legend(x= "topright" , legend = c("Alfa = 0.2", Q[, 1], "Alfa = 0.1", Q[, 2], "Alfa = 0.05", Q[, 3], "Alfa = 0.01", Q[, 4]))
    polygon(c(Den$x[Den$x <= value1], value1),
            c(Den$y[Den$x <= value1], 0),
            col = rgb(1, 0, 0, alpha = 0.2),
            border = 1)
    polygon(c(Den$x[Den$x <= value2], value2),
            c(Den$y[Den$x <= value2], 0),
            col = rgb(1, 0, 0, alpha = 0.2),
            border = 1)
    polygon(c(Den$x[Den$x <= value3], value3),
            c(Den$y[Den$x <= value3], 0),
            col = rgb(1, 0, 0, alpha = 0.2),
            border = 1)
    polygon(c(Den$x[Den$x <= value4], value4),
            c(Den$y[Den$x <= value4], 0),
            col = rgb(1, 0, 0, alpha = 0.2),
            border = 1)
    abline(v = mean(Tau_T), col = "red", lwd = 3, lty = 2)
  })
  
  
  output$Tau = renderPlot(
    Tau()
  )
  output$Taumu = renderPlot(
    Taumu()
    
  )
  output$TauTau = renderPlot(
    TauTau()
  )
}

shinyApp(ui=ui, server=server)