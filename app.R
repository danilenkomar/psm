Sys.setlocale("LC_ALL", "Russian")
library(shiny)
library(shinyjs)
library(bslib)
library(bsicons)
library(psych)
library(MatchIt)
library(dplyr)
library(purrr)
library(ggplot2)
library(survival)
# Определяем пользовательский интерфейс
ui <- fluidPage(
  
  # Заголовок приложения
  titlePanel("Propensity score matching"),
  
  # Вкладки с различными страницами
  tabsetPanel(
    id = "tabs",
    tabPanel("Загрузка данных", value = "tabMain",
             fluidRow(
               # Первая колонка с полем для загрузки CSV файла
               column(6, fileInput("file", "Выберите CSV файл"), 
                      textInput("covariates", "Введите названия ковариат (через запятую)"),
                      textInput("targets", "Введите названия целевых показателей (через запятую)"),
                      textInput("groups", "Введите название признака, определяющего группу для сравнения")),
               # Вторая колонка с блоком информации о датафрейме
               column(6, verbatimTextOutput("data_info"))
             ),
             # Поле для указания количества строк для отображения
             numericInput("n_rows", "Количество строк для отображения", value = 5),
             # Таблица для отображения первых n строк датафрейма
             card(tableOutput("data_table"))
    ),
    tabPanel("Ковариаты",id = "tabCov",   
             fluidRow(
               column(6,selectInput("methodSop", "Выберите метод сопоставления:",
                                                      c("Метод наилучшего соотвествия" = "nearest",
                                                        #"Метод генетического соответствия" = "genetic",
                                                        "Метод оптимального сопоставления" = "optimal"), selected = "optimal"),
             
             sliderInput("obs", "Для метода наилучшего соотвествия укажите максимальное расстояние (ст. откл.):",
                         min = 0.05, max = 1, value = 0.2, step = 0.05
             ),),
             column(3,verbatimTextOutput("data_info_psm"))
             ),
             fluidRow(
               # Первая колонка со сводной таблицей до PSM
               column(6,h4("Полный набор данных:"), card(tableOutput("cov_table"))),
               # Вторая колонка со сводной таблицей после PSM
               column(6, h4("После PSM:"), card(tableOutput("cov_table_after_PSM")))
             )),
    tabPanel("Целевые показатели", id = "tabTar",              
             fluidRow(
      # Первая колонка со сводной таблицей до PSM
      column(6,h4("Полный набор данных:"), card(tableOutput("tar_table"))),
      # Вторая колонка со сводной таблицей после PSM
      column(6, h4("После PSM:"), card(tableOutput("tar_table_after_PSM")))
    )),
    tabPanel("Визуализация",id = "tabVis",
             fluidRow(
               # Первая колонка с полем для загрузки CSV файла
               column(4, h3("Круговая диаграмма:"), 
                      selectInput("selected_bin", "Выберите бинарный признак:", choices = NULL)),
               # Вторая колонка с блоком информации о датафрейме
               column(4, h3("Гистограмма распределения:"), 
                      selectInput("selected_num_1", "Выберите числовой признак:", choices = NULL)
                      ),
               column(4, h3("Диаграмма размаха:"), 
                      selectInput("selected_num_2", "Выберите числовой признак:", choices = NULL),)
             ),
             h4("После PSM:"),
             fluidRow(
               # Первая колонка с полем для загрузки CSV файла
               column(4,
                      plotOutput("round_psm")),
               # Вторая колонка с блоком информации о датафрейме
               column(4, 
                      plotOutput("hist_psm")),
               column(4,
                      plotOutput("box_psm"))
             ),
             h4("Полный набор данных:"),
             fluidRow(
               # Первая колонка с полем для загрузки CSV файла
               column(4, 
                      plotOutput("round")),
               # Вторая колонка с блоком информации о датафрейме
               column(4, 
                      plotOutput("hist")),
               column(4, 
                      plotOutput("box"))
             )
    ),
    tabPanel("Оценка выживаемости",id = "tabSer", 
             textInput("year_oper", "Введите название признака, определяющего год проведения лечения"),
             textInput("year_die", "Введите название признака, определяющего год смерти"),
             textInput("cur_year", "Введите последний год сбора данных"),
             fluidRow(
               # Первая колонка со сводной таблицей до PSM
               column(6,h3("Полный набор данных:"),  plotOutput("sur")),
               # Вторая колонка со сводной таблицей после PSM
               column(6,h3("После PSM:"),  plotOutput("sur_psm"))
             ))
  )
)

# Определяем серверную часть
server <- function(input, output, session) {
  #Объявление функций---------  
  convert_to_numeric <- function(df) {
    for (col in names(df)) {
      if (is.numeric(df[[col]])) {
        unique_vals <- unique(df[[col]])
        if (length(unique_vals) > 2 || !all(unique_vals %in% c(0, 1))) {
          df[[col]] <- as.numeric(as.character(df[[col]]))
        }
      }
    }
    return(df)
  }
  get_numeric_columns <- function(data_frame) {
    numeric_columns <- sapply(data_frame, function(col) {
      is_numeric <- is.numeric(col)
      is_integer <- is.integer(col)
      if (is_numeric && !is_integer) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    })
    return(names(data_frame)[numeric_columns])
  }
  get_integer_columns <- function(data_frame) {
    integer_columns <- sapply(data_frame, function(col) {
      is_integer <- is.integer(col)
      if (is_integer) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    })
    return(names(data_frame)[integer_columns])
  }
  # Функция для вычисления абсолютных стандартизированных разностей
  calc_abs_std_diff <- function(mgp, bgp) {
    round(((mean(mgp, na.rm = TRUE) - mean(bgp, na.rm = TRUE)) / sqrt((var(mgp, na.rm = TRUE) + var(bgp, na.rm = TRUE)) / 2)), 4)
  }
  
  # Функция для вычисления медианы, квартилей и процента встречаемости признака
  summary_stats <- function(data) {
    med <- median(data, na.rm = TRUE)
    q25 <- quantile(data, 0.25,  na.rm = TRUE)
    q75 <- quantile(data, 0.75,  na.rm = TRUE)
    
    if (class(data) == "integer") {
      unique_perc <- round((sum(data == 1,  na.rm = TRUE) / sum(!is.na(data))) * 100, 2)
      result <- paste0(unique_perc, "%")
    } else {
      result <- paste0(med,' (', q25,';', q75,')')
    }
    
    return(result)
  }
  
  #Функция получения сводной таблицы
  sum_tab_get <-function(mgp, bgp, feature){
    mgp_f <- mgp[,feature]
    bgp_f <- bgp[,feature]
    # Создание таблицы с характеристиками для первого набора данных
    summary_table_mgp <- t(summarise_all(mgp_f, summary_stats))
    colnames(summary_table_mgp) <- "Группа 1"
    
    # Создание таблицы с характеристиками для второго набора данных
    summary_table_bgp <- t(summarise_all(bgp_f, summary_stats))
    colnames(summary_table_bgp) <- "Группа 0"
    # Вычисление абсолютных стандартизированных разностей
    abs_std_diff <- map2_dbl(mgp_f, bgp_f, calc_abs_std_diff)
    
    # Добавление абсолютных стандартизированных разностей в таблицу
    summary_table <- cbind(summary_table_mgp,summary_table_bgp, abs_std_diff)
    
    # Вычисление p-values с помощью критерия Манна-Уитни для количественных признаков и критерия Фишера для бинарных
    p <- map2_dbl(mgp_f, bgp_f, function(x, y) {
      if (class(data) == "integer") {
        round(fisher.test(table(x, y))$p.value,4)
      } else {
        round(wilcox.test(x, y)$p.value,4)
        #round(t.test(x, y)$p.value,4)
      }
    })
    summary_table <- cbind(summary_table, p)
    colnames(summary_table)[colnames(summary_table) == "abs_std_diff"] <- "Абсолютные_ст._разности"
    return(summary_table)
  }
  #Функция построения кривой выживаемости
  surv<- function(op, group,time_mun,time_dead,current){
    time <- op[[time_dead]] - op[[time_mun]]
    event <- !is.na(op[[time_dead]]) # Пациент умер, если время смерти известно
    
    # Учитываем случаи, когда время смерти неизвестно, но пациенты еще живы
    # Пусть время выживания для таких пациентов будет текущему времени минус год операции
    time[is.na(op[[time_dead]])] <- current - op[is.na(op[[time_dead]]),time_mun ] 
    event[is.na(op[[time_dead]])] <- FALSE
    print(event)
    surv_obj <- Surv(time = time, event = event)
    formula <- as.formula(paste("surv_obj ~", group))
    # Построение кривой выживаемости методом Каплана-Майера для обеих групп
    km_curve <- survfit(formula, data = op)
    # Проведение лог-рангового теста для сравнения двух групп
    logrank_test <- survdiff(formula =formula, data = op)
    p_value <- 1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1)
    print(p_value)
    # Визуализация кривой выживаемости
    par(mfrow=c(1, 1))
    plot(km_curve, xlab = "Год", ylab = "Доля выживших", col = c("blue", "red"), lwd = 2)
    legend("bottomleft", legend = c("Группа 0", "Группа 1"), col = c("blue", "red"), lwd = 2)
    text(x = max(km_curve$time), y = 0, labels = paste("p = ", round(p_value,4) ), adj = c(1, 0))
  }
  
  #Гистограмма
  hist_view <- function(mgp, bgp, col_name){
    par(mfrow = c(1, 1))
    hist(mgp[[col_name]], breaks=25, col=rgb(0,0,0,alpha=.5), xlab="Graph", main=col_name,probability = TRUE) #нормировка
    hist(bgp[[col_name]], breaks=25, col=rgb(0,0,1,alpha=.5), xlab="Graph", add = TRUE,probability = TRUE) 
    legend("topright", c("Группа 1","Группа 0"), cex = 0.8, fill =c(rgb(0,0,0,alpha=.5),rgb(0,0,1,alpha=.5)))
    
  }
  #Диаграмма с усиками 
  box_view <- function(data, group, col_name){
    par(mfrow = c(1, 1))
    ggplot(data, aes_string(x = paste0("factor(", group, ")"), y = col_name, fill = paste0("factor(", group, ")"))) +
      geom_boxplot() +
      theme_minimal() +
      labs(title = col_name, x = "Группа", y = "Значение")
  }
  #Круговая диаграмма
  pie_view <- function(mgp, bgp, col_name){
    par(mfrow=c(1, 2))
    count_0 <- sum(mgp[col_name] == 0)
    count_1 <- sum(mgp[col_name] == 1)
    print (count_0)
    print (nrow(mgp))
    # Преобразование в проценты
    percent_0 <- round(((count_0 / nrow(mgp)) * 100),1)
    percent_1 <- round(((count_1 / nrow(mgp)) * 100),1)
    
    x<-c(count_0,count_1)
    piepercent<-c(percent_0, percent_1)
    pie(x, piepercent, radius=1, main="Группа 1", col=c("blue", "red"), clockwise=TRUE)
    legend("topright", c(paste0(col_name," = 0"), paste0(col_name," = 1")), cex = 0.8, fill =c("blue", "red"))
    
    count_0 <- sum(bgp[col_name] == 0)
    count_1 <- sum(bgp[col_name] == 1)
    
    # Преобразование в проценты
    percent_0 <- round(((count_0 / nrow(bgp)) * 100),1)
    percent_1 <- round(((count_1 / nrow(bgp)) * 100),1)
    
    x<-c(count_0,count_1)
    piepercent<-c(percent_0, percent_1)
    pie(x, piepercent, radius=1, main="Группа 0", col=c("blue", "red"), clockwise=TRUE)
    legend("topright", c(paste0(col_name," = 0"), paste0(col_name," = 1")), cex = 0.8, fill =c("blue", "red"))
    
  }
  #Серверная часть обработки-----
  df1<- reactiveVal(NULL)
  df1_PSM<- reactiveVal(NULL)
  cov<- reactiveVal(NULL)
  tar<- reactiveVal(NULL)
  gro<- reactiveVal(NULL)
  bin<- reactiveVal(NULL)
  num<- reactiveVal(NULL)
  year_oper<- reactiveVal(NULL)
  year_die<- reactiveVal(NULL)
  year_cur<- reactiveVal(NULL)
  
  observeEvent(input$tabs, {
    print (is.null(df1()))
    if (is.null(df1())||is.null(cov())||is.null(tar())||is.null(gro())) {
      updateTabsetPanel(session = getDefaultReactiveDomain(), "tabs",
                        selected = "tabMain")
    }
    else{
      formula <- as.formula(paste(gro()," ~", paste(cov(), collapse = " + ")))
      if (input$methodSop == "optimal"){
        matched_pairs <- matchit(formula, data = df1(),
                                     method = "optimal", distance = "glm", link = "logit")
        
        df_PSM  <- match.data(matched_pairs)
        df1_PSM(df_PSM)
        output$cov_table_after_PSM<-renderTable({sum_tab_get(df_PSM[df_PSM[[gro()]] == 1,], df_PSM[df_PSM[[gro()]] == 0,], cov())}, rownames = TRUE)
      }
      else{
        matched_pairs <- matchit(formula, data = df1(),
                                 method = input$methodSop, caliper = input$obs, distance = "glm")
        
        df_PSM  <- match.data(matched_pairs)
        df1_PSM(df_PSM)
        output$cov_table_after_PSM<-renderTable({sum_tab_get(df_PSM[df_PSM[[gro()]] == 1,], df_PSM[df_PSM[[gro()]] == 0,], cov())}, rownames = TRUE)
      }
      output$tar_table_after_PSM<-renderTable({sum_tab_get(df_PSM[df_PSM[[gro()]] == 1,], df_PSM[df_PSM[[gro()]] == 0,], tar())}, rownames = TRUE)
      df = df1()
      output$data_info_psm <- renderPrint({
        cat("Количество строк в группе 1 до PSM:", nrow(df[df[[gro()]] == 1,]), "\n")
        cat("Количество строк в группе 0 до PSM:", nrow(df[df[[gro()]] == 0,]), "\n")
        cat("Количество строк в группе 1 после PSM:", nrow(df_PSM[df_PSM[[gro()]] == 1,]), "\n")
        cat("Количество строк в группе 0 после PSM:", nrow(df_PSM[df_PSM[[gro()]] == 0,]), "\n")
      })
      output$cov_table<-renderTable({sum_tab_get(df[df[[gro()]] == 1,], df[df[[gro()]] == 0,], cov())}, rownames = TRUE)
      output$tar_table<-renderTable({sum_tab_get(df[df[[gro()]] == 1,], df[df[[gro()]] == 0,], tar())}, rownames = TRUE)
    }
  })
  # Функция для чтения загруженного CSV файла и вывода первых n строк

  observeEvent(input$file, {
    req(input$file) # Требуем наличие загруженного файла
    df <- read.csv(input$file$datapath, encoding = "UTF-8")
    df <- convert_to_numeric(df)
    bin(get_integer_columns(df))
    num(get_numeric_columns(df))
    df1(df)
    output$data_info <- renderPrint({
      cat("Количество строк в датафрейме:", nrow(df), "\n")
      cat("Количество признаков в датафрейме:", ncol(df), "\n")
      cat("Признаки:", "\n")
      col_names <- names(df)
      num_cols <- length(col_names)
      cols_per_block <- 3
      num_blocks <- ceiling(num_cols / cols_per_block)
      for (i in 1:num_blocks) {
        start_index <- (i - 1) * cols_per_block + 1
        end_index <- min(i * cols_per_block, num_cols)
        cat(paste(col_names[start_index:end_index], "-", sapply(df[start_index:end_index], class), collapse = "; "), "\n")
      }
    })
    updateTextInput(session, "n_rows", value = min(5, nrow(df)))
    output$data_table <- renderTable({
      head(df, n = input$n_rows) # Выводим первые n строк датафрейма
    },rownames = TRUE)
    updateSelectInput(session, "selected_bin", choices = bin()) #Оставить в наборе только признаки, указанные как целевые или ковариаты??
    updateSelectInput(session, "selected_num_1", choices = num())
    updateSelectInput(session, "selected_num_2", choices = num())
    
  })
  
  observeEvent(input$targets, {
    req(input$targets)  # Проверяем, что значение введено
    tar(unlist(strsplit(input$targets, ", "))) #сразу разбиваем список
  })
  
  observeEvent(input$groups, {
    req(input$groups)  # Проверяем, что значение введено
    gro(input$groups)
  })
  
  observeEvent(input$year_oper , {
    req(input$year_oper)  # Проверяем, что значение введено
    year_oper(input$year_oper)
  })
  
  observeEvent(input$year_die, {
    req(input$year_die)  # Проверяем, что значение введено
    year_die(input$year_die)
  })
  
  observeEvent(input$cur_year, {
    req(input$cur_year)  # Проверяем, что значение введено
    year_cur(input$cur_year)
  })
  
  observeEvent(input$covariates, {
    req(input$covariates)  # Проверяем, что значение введено
    cov(unlist(strsplit(input$covariates, ", "))) #сразу разбиваем список
  })
  
  observeEvent(input$obs, {
    req(input$obs)  # Проверяем, что значение введено
    if (!(is.null(df1())||is.null(cov())||is.null(tar())||is.null(gro()))){
    formula <- as.formula(paste(gro()," ~", paste(cov(), collapse = " + ")))
    if (input$methodSop == "optimal"){
      matched_pairs <- matchit(formula, data = df1(),
                               method = "optimal", distance = "glm", link = "logit")
      
      df_PSM  <- match.data(matched_pairs)
      df1_PSM(df_PSM)
      output$cov_table_after_PSM<-renderTable({sum_tab_get(df_PSM[df_PSM[[gro()]] == 1,], df_PSM[df_PSM[[gro()]] == 0,], cov())}, rownames = TRUE)
    }
    else{
      matched_pairs <- matchit(formula, data = df1(),
                               method = input$methodSop, caliper = input$obs, distance = "glm")
      
      df_PSM  <- match.data(matched_pairs)
      df1_PSM(df_PSM)
      output$cov_table_after_PSM<-renderTable({sum_tab_get(df_PSM[df_PSM[[gro()]] == 1,], df_PSM[df_PSM[[gro()]] == 0,], cov())}, rownames = TRUE)
    }
    output$tar_table_after_PSM<-renderTable({sum_tab_get(df_PSM[df_PSM[[gro()]] == 1,], df_PSM[df_PSM[[gro()]] == 0,], tar())}, rownames = TRUE)
    df = df1()
    output$data_info_psm <- renderPrint({
      cat("Количество строк в группе 1 до PSM:", nrow(df[df[[gro()]] == 1,]), "\n")
      cat("Количество строк в группе 0 до PSM:", nrow(df[df[[gro()]] == 0,]), "\n")
      cat("Количество строк в группе 1 после PSM:", nrow(df_PSM[df_PSM[[gro()]] == 1,]), "\n")
      cat("Количество строк в группе 0 после PSM:", nrow(df_PSM[df_PSM[[gro()]] == 0,]), "\n")
    })}
  })
  
  observeEvent(input$methodSop, {
    req(input$methodSop)  # Проверяем, что значение введено
    if (!(is.null(df1())||is.null(cov())||is.null(tar())||is.null(gro()))){
    formula <- as.formula(paste(gro()," ~", paste(cov(), collapse = " + ")))
    if (input$methodSop == "optimal"){
      matched_pairs <- matchit(formula, data = df1(),
                               method = "optimal", distance = "glm", link = "logit")
      
      df_PSM  <- match.data(matched_pairs)
      df1_PSM(df_PSM)
      output$cov_table_after_PSM<-renderTable({sum_tab_get(df_PSM[df_PSM[[gro()]] == 1,], df_PSM[df_PSM[[gro()]] == 0,], cov())}, rownames = TRUE)
    }
    else{
      matched_pairs <- matchit(formula, data = df1(),
                               method = input$methodSop, caliper = input$obs, distance = "glm")
      
      df_PSM  <- match.data(matched_pairs)
      df1_PSM(df_PSM)
      output$cov_table_after_PSM<-renderTable({sum_tab_get(df_PSM[df_PSM[[gro()]] == 1,], df_PSM[df_PSM[[gro()]] == 0,], cov())}, rownames = TRUE)
    }
    output$tar_table_after_PSM<-renderTable({sum_tab_get(df_PSM[df_PSM[[gro()]] == 1,], df_PSM[df_PSM[[gro()]] == 0,], tar())}, rownames = TRUE)
    df = df1()
    output$data_info_psm <- renderPrint({
      cat("Количество строк в группе 1 до PSM:", nrow(df[df[[gro()]] == 1,]), "\n")
      cat("Количество строк в группе 0 до PSM:", nrow(df[df[[gro()]] == 0,]), "\n")
      cat("Количество строк в группе 1 после PSM:", nrow(df_PSM[df_PSM[[gro()]] == 1,]), "\n")
      cat("Количество строк в группе 0 после PSM:", nrow(df_PSM[df_PSM[[gro()]] == 0,]), "\n")
    })}
  })
  observeEvent(input$selected_bin, {
    print(input$selected_bin)
    if (input$selected_bin != ""){
    
    output$round <- renderPlot({
      df = df1()
      pie_view(df[df[[gro()]] == 1,],df[df[[gro()]] == 0,],input$selected_bin)})
    output$round_psm <- renderPlot({
      df_psm = df1_PSM()
      pie_view(df_psm[df_psm[[gro()]] == 1,],df_psm[df_psm[[gro()]] == 0,],input$selected_bin)})}
  })
  observeEvent(input$selected_num_1, {
    print(input$selected_num_1)
    if (input$selected_num_1 != ""){
      
      output$hist <- renderPlot({
        df = df1()
        hist_view(df[df[[gro()]] == 1,],df[df[[gro()]] == 0,],input$selected_num_1)})
      output$hist_psm <- renderPlot({
        df_psm = df1_PSM()
        hist_view(df_psm[df_psm[[gro()]] == 1,],df_psm[df_psm[[gro()]] == 0,],input$selected_num_1)})}
  })
  observeEvent(input$selected_num_2, {
    print(input$selected_num_2)
    if (input$selected_num_2 != ""){
      
      output$box <- renderPlot({
        df = df1()
        box_view(df, gro() ,input$selected_num_2)})
      output$box_psm <- renderPlot({
        df_psm = df1_PSM()
        box_view(df_psm, gro() ,input$selected_num_2)})}
  })
  output$sur <- renderPlot({
    if(!is.null(year_oper())&&!is.null(year_die())&&!is.null(year_cur())){
    df = df1()
    surv(df, gro(),input$year_oper,input$year_die,as.numeric(input$cur_year))
    }})
  
  output$sur_psm <- renderPlot({
    if(!is.null(year_oper())&&!is.null(year_die())&&!is.null(year_cur())){
      df = df1_PSM()
      surv(df, gro(),input$year_oper,input$year_die,as.numeric(input$cur_year))
    }})
  
}
# Запускаем приложение Shiny
shinyApp(ui = ui, server = server)