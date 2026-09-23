library(data.table)

files <- c(
"SRR3142259.count",
"SRR3142260.count",
"SRR3142266.count",
"SRR3142261.count",
"SRR3142267.count",
"SRR3142262.count",
"SRR3142268.count",
"SRR3142263.count",
"SRR3142265.count",
"SRR3142269.count"
)

############################
# 读取第一个文件
############################

x <- fread(files[1], skip=1)

counts <- x[, .(
    Geneid,
    Length,
    SRR3142259 = get(names(x)[7])
)]

############################
# 合并其余样本
############################

for(i in 2:length(files)){

    f <- files[i]

    sample <- sub(".count","",f)

    x <- fread(f, skip=1)

    tmp <- x[, .(
        Geneid,
        Count = get(names(x)[7])
    )]

    setnames(tmp,"Count",sample)

    counts <- merge(
        counts,
        tmp,
        by="Geneid"
    )
}

############################
# TPM计算
############################

count.mat <- as.matrix(
    counts[, -(1:2)]
)

length.kb <- counts$Length / 1000

rpk <- sweep(
    count.mat,
    1,
    length.kb,
    "/"
)

scale.factor <- colSums(rpk) / 1e6

tpm <- sweep(
    rpk,
    2,
    scale.factor,
    "/"
)

tpm.df <- data.frame(
    Geneid = counts$Geneid,
    tpm,
    check.names = FALSE
)

############################
# 输出
############################

write.table(
    tpm.df,
    file="merged_TPM.txt",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("TPM calculation finished!\n")
