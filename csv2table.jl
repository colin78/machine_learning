using Mustache
using DataFrames

format_cell(x) = "$x"
format_cell(x::String) = replace(x, "%", "\\%")
format_cell(x::NAtype) = "-"
format_cell(x::Integer) = "$x"
# format_cell(x::Integer) = x == 0 ? @sprintf("%.4f", x) : "$x"
format_cell(x::Real) = @sprintf("%.4f", x)

function df_to_table(df::DataFrame, format, header, label, caption,
                     position="tb")
  row = join(["{{$x}}" for x in map(string, names(df))], " & ")

  tpl = """
        \\begin{table}[$position]
        \\centering
        \\begin{tabular}{$format}
          \\toprule
          $header\\\\
          \\midrule
        {{#df}}  $row \\\\
        {{/df}}  \\bottomrule
        \\end{tabular}
        \\caption{$caption}\\label{tab:$label}
        \\end{table}
        """

  render(tpl, ["df"=>df])
end

# Example: maketable("results\\V","ccc","V1 & V2 & V3","label","caption3")
function maketable(filename, format, header, label, caption; position="tb", headerPresent=false)
  df = readtable("$filename.csv", header=headerPresent)
  # Convert all numbers to 4dp precision (except non-zero ints)
  for i in 1:length(df)
    df[i] = map(format_cell, df[i])
  end

  f = open("$filename.tex", "w")
  write(f, df_to_table(df, format, header, label, caption, position))
  close(f)
end

function makecomparisontable(filename, format, header, label, caption)
  df = readtable("$filename.csv", header=true)
  n = length(df[1])
  numresults = length(df)

  inds_lr   = Array(Any, n)
  inds_dsvm = Array(Any, n)
  inds_cart = Array(Any, n)
  inds_all  = Array(Any, n)

  # Convert all numbers to 4dp precision (except non-zero ints)
  for i in 1:n
    max_lr   = max(df[i, 4], df[i, 5], df[i, 6], df[i, 7])
    max_dsvm = max(df[i, 8], df[i, 9], df[i, 10], df[i, 11])
    max_cart = max(df[i, 12], df[i, 13], df[i, 14], df[i, 15])
    max_all  = max(max_lr, max_dsvm, max_cart)
    inds_lr[i] = filter(x -> (isna(df[i, x]) ? -10 : df[i, x]) == max_lr,
                        4:7)
    inds_dsvm[i] = filter(x -> (isna(df[i, x]) ? -10 : df[i, x]) == max_dsvm,
                          8:11)
    inds_cart[i] = filter(x -> (isna(df[i, x]) ? -10 : df[i, x]) == max_cart,
                          12:15)
    inds_all[i] = filter(x -> (isna(df[i, x]) ? -10 : df[i, x]) == max_all,
                         4:15)
  end

  for j in 1:numresults
    df[j] = map(format_cell, df[j])
  end

  for i in 1:n
    for j in 1:numresults
      if j in inds_lr[i]
        df[i, j] = "\\textbf{$(df[i, j])}"
      end
      if j in inds_dsvm[i]
        df[i, j] = "\\textbf{$(df[i, j])}"
      end
      if j in inds_cart[i]
        df[i, j] = "\\textbf{$(df[i, j])}"
      end
      if j in inds_all[i]
        df[i, j] = "\\underline{$(df[i, j])}"
      end
    end
  end

  df = df[[1:3; 8:11; 4:7; 12:15]]

  row = join(["{{$x}}" for x in map(string, names(df))], " & ")

  tpl = """
        {\\tiny
        \\begin{longtable}{$format}
          \\toprule
          $header\\\\
          \\midrule
          \\endfirsthead
          \\toprule
          $header\\\\
          \\midrule
          \\endhead
          \\bottomrule
          \\\\ \\caption{$caption}\\label{tab:$label}\\vspace{1mm}
          \\endfirstfoot
          \\bottomrule
          \\\\ \\caption{(Cont.) $caption}\\vspace{1mm}
          \\endfoot
          {{#df}}  $row\\\\
          {{/df}}
        \\end{longtable}}
        """

  data = render(tpl, ["df"=>df])

  f = open("$filename.tex", "w")
  write(f, data)
  close(f)
end

function makecomparisontable_state(filename, format, header, label, caption;
                                   small=false)
  df = readtable("$filename.csv", header=true)
  n = length(df[1])
  numresults = length(df)

  # data = Matrix(df[:, 4:end]) # assumes everything is float
  inds_lr   = Array(Any, n)
  inds_dsvm = Array(Any, n)
  inds_cart = Array(Any, n)
  inds_all  = Array(Any, n)

  for i in 1:n
    max_lr   = max(df[i, 4], df[i, 5])
    max_dsvm = max(df[i, 6], df[i, 7])
    max_cart = max(df[i, 8], df[i, 9])
    max_all  = max(max_lr, max_dsvm, max_cart)
    inds_lr[i] = filter(x -> (isna(df[i, x]) ? -10 : df[i, x]) == max_lr,
                        4:5)
    inds_dsvm[i] = filter(x -> (isna(df[i, x]) ? -10 : df[i, x]) == max_dsvm,
                          6:7)
    inds_cart[i] = filter(x -> (isna(df[i, x]) ? -10 : df[i, x]) == max_cart,
                          8:9)
    inds_all[i] = filter(x -> (isna(df[i, x]) ? -10 : df[i, x]) == max_all,
                         4:9)
  end

  # Convert all numbers to 4dp precision (except non-zero ints)
  for j in 1:numresults
    df[j] = map(format_cell, df[j])
  end

  for i in 1:n
    for j in 1:numresults
      if j in inds_lr[i]
        df[i, j] = "\\textbf{$(df[i, j])}"
      end
      if j in inds_dsvm[i]
        df[i, j] = "\\textbf{$(df[i, j])}"
      end
      if j in inds_cart[i]
        df[i, j] = "\\textbf{$(df[i, j])}"
      end
      # if j in inds_all[i]
      #   df[i, j] = "\\underline{$(df[i, j])}"
      # end
    end
  end

  df = df[[1:3; 6:7; 4:5; 8:9]]
  if small
    df = df[1:7]
    filename = string(filename, "_small")
  end

  row = join(["{{$x}}" for x in map(string, names(df))], " & ")

  tpl = """
        {\\scriptsize
        \\begin{longtable}{$format}
          \\toprule
          $header\\\\
          \\midrule
          \\endfirsthead
          \\toprule
          $header\\\\
          \\midrule
          \\endhead
          \\bottomrule
          \\\\ \\caption{$caption}\\label{tab:$label}\\vspace{1mm}
          \\endfirstfoot
          \\bottomrule
          \\\\ \\caption{(Cont.) $caption}\\vspace{1mm}
          \\endfoot
          {{#df}}  $row\\\\
          {{/df}}
        \\end{longtable}}
        """

  data = render(tpl, ["df"=>df])

  f = open("$filename.tex", "w")
  write(f, data)
  close(f)
end


function makewincounttable(filename, label, caption)
  df = readtable("$filename.csv", header=false)
  table = """
          \\begin{table}[tb]
          \\centering
          \\begin{tabular}{ccccc}
            \\toprule
            Nominal Method & Robustness Type & Wins & Losses & Ties\\\\
            \\midrule
            \\multirow{3}{*}{SVM} & Features & $(df[1,3]) & $(df[2,3]) & $(df[3,3])\\\\
            & Labels & $(df[4,3]) & $(df[5,3]) & $(df[6,3])\\\\
            & Both & $(df[7,3]) & $(df[8,3]) & $(df[9,3])\\\\
            \\midrule
            \\multirow{3}{*}{Logistic Regression} & Features & $(df[10,3]) & $(df[11,3]) & $(df[12,3])\\\\
            & Labels & $(df[13,3]) & $(df[14,3]) & $(df[15,3])\\\\
            & Both & $(df[16,3]) & $(df[17,3]) & $(df[18,3])\\\\
            \\midrule
            \\multirow{3}{*}{CART} & Features & $(df[19,3]) & $(df[20,3]) & $(df[21,3])\\\\
            & Labels & $(df[22,3]) & $(df[23,3]) & $(df[24,3])\\\\
            & Both & $(df[25,3]) & $(df[26,3]) & $(df[27,3])\\\\
            \\bottomrule
          \\end{tabular}
          \\caption{$caption}\\label{tab:$label}
          \\end{table}
          """
  f = open("$filename.tex", "w")
  write(f, table)
  close(f)
end

function makeacctable(prefix, label, caption)
  df_logit = readtable("$(prefix)_logit.csv", header=false)
  df_svm = readtable("$(prefix)_svm.csv", header=false)
  df_cart = readtable("$(prefix)_cart.csv", header=false)

  table = """
          \\begin{table}[tb]
          \\centering
          \\begin{tabular}{cccccc}
            \\toprule
            Nominal Method & Nominal Accuracy & Wins & Losses & Ties & Robust Improvement\\\\
            \\midrule
            \\multirow{5}{*}{SVM} & $(format_cell(df_svm[1,1])) & $(df_svm[1,2]) & $(df_svm[1,3]) & $(df_svm[1,4]) & $(format_cell(df_svm[1,5]))\\\\
            & $(format_cell(df_svm[2,1])) & $(df_svm[2,2]) & $(df_svm[2,3]) & $(df_svm[2,4]) & $(format_cell(df_svm[2,5]))\\\\
            & $(format_cell(df_svm[3,1])) & $(df_svm[3,2]) & $(df_svm[3,3]) & $(df_svm[3,4]) & $(format_cell(df_svm[3,5]))\\\\
            & $(format_cell(df_svm[4,1])) & $(df_svm[4,2]) & $(df_svm[4,3]) & $(df_svm[4,4]) & $(format_cell(df_svm[4,5]))\\\\
            & $(format_cell(df_svm[5,1])) & $(df_svm[5,2]) & $(df_svm[5,3]) & $(df_svm[5,4]) & $(format_cell(df_svm[5,5]))\\\\
            \\midrule
            \\multirow{5}{*}{Logistic Regression} & $(format_cell(df_logit[1,1])) & $(df_logit[1,2]) & $(df_logit[1,3]) & $(df_logit[1,4]) & $(format_cell(df_logit[1,5]))\\\\
            & $(format_cell(df_logit[2,1])) & $(df_logit[2,2]) & $(df_logit[2,3]) & $(df_logit[2,4]) & $(format_cell(df_logit[2,5]))\\\\
            & $(format_cell(df_logit[3,1])) & $(df_logit[3,2]) & $(df_logit[3,3]) & $(df_logit[3,4]) & $(format_cell(df_logit[3,5]))\\\\
            & $(format_cell(df_logit[4,1])) & $(df_logit[4,2]) & $(df_logit[4,3]) & $(df_logit[4,4]) & $(format_cell(df_logit[4,5]))\\\\
            & $(format_cell(df_logit[5,1])) & $(df_logit[5,2]) & $(df_logit[5,3]) & $(df_logit[5,4]) & $(format_cell(df_logit[5,5]))\\\\
            \\midrule
             \\multirow{5}{*}{CART} & $(format_cell(df_cart[1,1])) & $(df_cart[1,2]) & $(df_cart[1,3]) & $(df_cart[1,4]) & $(format_cell(df_cart[1,5]))\\\\
            & $(format_cell(df_cart[2,1])) & $(df_cart[2,2]) & $(df_cart[2,3]) & $(df_cart[2,4]) & $(format_cell(df_cart[2,5]))\\\\
            & $(format_cell(df_cart[3,1])) & $(df_cart[3,2]) & $(df_cart[3,3]) & $(df_cart[3,4]) & $(format_cell(df_cart[3,5]))\\\\
            & $(format_cell(df_cart[4,1])) & $(df_cart[4,2]) & $(df_cart[4,3]) & $(df_cart[4,4]) & $(format_cell(df_cart[4,5]))\\\\
            & $(format_cell(df_cart[5,1])) & $(df_cart[5,2]) & $(df_cart[5,3]) & $(df_cart[5,4]) & $(format_cell(df_cart[5,5]))\\\\
            \\bottomrule
          \\end{tabular}
          \\caption{$caption}\\label{tab:$label}
          \\end{table}
          """
  f = open("$prefix.tex", "w")
  write(f, table)
  close(f)
end



function makeacctable_by_nomacc_p(prefix, label, caption)
  df_logit = readtable("$(prefix)_logit.csv", header=false)
  df_svm = readtable("$(prefix)_svm.csv", header=false)
  df_cart_mip = readtable("$(prefix)_cart_mip.csv", header=false)
  df_cart = readtable("$(prefix)_cart.csv", header=false)

  table = """
          \\begin{table}[tb]
          \\centering
          \\begin{tabular}{cccccc}
            \\toprule
            Baseline Method & Region & Wins & Losses & Ties & Robust Improvement\\\\
            \\midrule
            \\multirow{2}{*}{\\parbox{5cm}{\\centering Nominal SVM}} & $(format_cell(df_svm[1,1])) & $(df_svm[1,2]) & $(df_svm[1,3]) & $(df_svm[1,4]) & $(format_cell(df_svm[1,5]))\\\\
            & $(format_cell(df_svm[2,1])) & $(df_svm[2,2]) & $(df_svm[2,3]) & $(df_svm[2,4]) & $(format_cell(df_svm[2,5]))\\\\
            \\midrule
            \\multirow{2}{*}{\\parbox{5cm}{\\centering Nominal\\\\ Logistic Regression}} & $(format_cell(df_logit[1,1])) & $(df_logit[1,2]) & $(df_logit[1,3]) & $(df_logit[1,4]) & $(format_cell(df_logit[1,5]))\\\\
            & $(format_cell(df_logit[2,1])) & $(df_logit[2,2]) & $(df_logit[2,3]) & $(df_logit[2,4]) & $(format_cell(df_logit[2,5]))\\\\
            \\midrule
            \\multirow{2}{*}{\\parbox{5cm}{\\centering Nominal\\\\ Optimal Decision Trees}} & $(format_cell(df_cart_mip[1,1])) & $(df_cart_mip[1,2]) & $(df_cart_mip[1,3]) & $(df_cart_mip[1,4]) & $(format_cell(df_cart_mip[1,5]))\\\\
            & $(format_cell(df_cart_mip[2,1])) & $(df_cart_mip[2,2]) & $(df_cart_mip[2,3]) & $(df_cart_mip[2,4]) & $(format_cell(df_cart_mip[2,5]))\\\\
            \\midrule
            \\multirow{2}{*}{CART} & $(format_cell(df_cart[1,1])) & $(df_cart[1,2]) & $(df_cart[1,3]) & $(df_cart[1,4]) & $(format_cell(df_cart[1,5]))\\\\
            & $(format_cell(df_cart[2,1])) & $(df_cart[2,2]) & $(df_cart[2,3]) & $(df_cart[2,4]) & $(format_cell(df_cart[2,5]))\\\\
            \\bottomrule
          \\end{tabular}
          \\caption{$caption}\\label{tab:$label}
          \\end{table}
          """
  f = open("$prefix.tex", "w")
  write(f, table)
  close(f)
end



function makeacctable_by_nomacc_p_SOA(prefix, label, caption)
  df_logit_reg = readtable("$(prefix)_logit_reg.csv", header=false)
  df_svm_classical = readtable("$(prefix)_svm_sk.csv", header=false)
  df_cart = readtable("$(prefix)_cart.csv", header=false)

  table = """
          \\begin{table}[tb]
          \\centering
          \\begin{tabular}{cccccc}
            \\toprule
            Baseline Method & Region & Wins & Losses & Ties & Robust Improvement\\\\
            \\midrule
            \\multirow{2}{*}{Classical SVM} & $(format_cell(df_svm_classical[1,1])) & $(df_svm_classical[1,2]) & $(df_svm_classical[1,3]) & $(df_svm_classical[1,4]) & $(format_cell(df_svm_classical[1,5]))\\\\
            & $(format_cell(df_svm_classical[2,1])) & $(df_svm_classical[2,2]) & $(df_svm_classical[2,3]) & $(df_svm_classical[2,4]) & $(format_cell(df_svm_classical[2,5]))\\\\
            \\midrule
            \\multirow{2}{*}{\\parbox{5cm}{\\centering Regularized\\\\ Logistic Regression}} & $(format_cell(df_logit_reg[1,1])) & $(df_logit_reg[1,2]) & $(df_logit_reg[1,3]) & $(df_logit_reg[1,4]) & $(format_cell(df_logit_reg[1,5]))\\\\
            & $(format_cell(df_logit_reg[2,1])) & $(df_logit_reg[2,2]) & $(df_logit_reg[2,3]) & $(df_logit_reg[2,4]) & $(format_cell(df_logit_reg[2,5]))\\\\
            \\bottomrule
          \\end{tabular}
          \\caption{$caption}\\label{tab:$label}
          \\end{table}
          """
  f = open("$(prefix)_SOA.tex", "w")
  write(f, table)
  close(f)
end
