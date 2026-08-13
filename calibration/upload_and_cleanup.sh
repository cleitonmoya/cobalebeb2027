#!/bin/bash
# upload_and_cleanup.sh
#
# Sobe uma lista de arquivos .rds pro OneDrive (via rclone), em lotes
# paralelos, e SO apaga cada arquivo do cluster depois de confirmar (via
# "rclone check", comparando checksum local x remoto) que a copia chegou
# integra -- nunca apaga so porque "rclone copy nao deu erro".
#
# Retomavel: se um arquivo ja existir no remoto com checksum batendo
# (ex.: rodada anterior que ja tinha subido, ou o script foi interrompido
# no meio), pula o upload e vai direto pra verificacao/remocao -- seguro
# rodar de novo do zero sem duplicar trabalho.
#
# Log em tempo real (uma linha por evento, escrita atomica -- linhas
# curtas via >> sao seguras mesmo com varios processos em paralelo
# escrevendo no mesmo arquivo, sem precisar de flock).
#
# Uso:
#   bash upload_and_cleanup.sh
# (edite REMOTE_DIR e a lista FILES abaixo antes de rodar)

set -u
shopt -s nullglob

REMOTE_DIR="onedrive:calibracao_backup/"
LOG="$HOME/upload_and_cleanup.log"
PARALLEL=8

# Descoberta DINAMICA dos arquivos, em vez de uma lista fixa transcrita a
# mao -- mais seguro. (Uma tentativa de reconciliar manualmente a lista
# completa contra o que ja estava no OneDrive bateu 63 arquivos em vez de
# 64, e um arquivo apareceu como "ja enviado" sem constar na listagem
# local copiada -- risco real de erro de transcricao dos dois lados.) A
# funcao upload_verify_delete() ja verifica, PARA CADA ARQUIVO, se ele
# already existe no remoto com checksum batendo antes de decidir subir --
# entao rodar sobre TODOS os .rds do diretorio atual (nao so os
# "faltantes" pre-calculados) e seguro e correto: os que ja foram
# enviados corretamente sao pulados automaticamente (ver o log "JA ESTAVA
# NO REMOTO"), sem re-upload nem duplicacao de trabalho.
FILES=(*.rds)

log() {
    # Data/hora + mensagem, uma linha, flush imediato (>> ja e sem buffer
    # -- cada chamada abre/escreve/fecha o arquivo).
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >> "$LOG"
}

upload_verify_delete() {
    local f="$1"

    if [ ! -f "$f" ]; then
        log "PULADO (nao existe localmente, provavelmente ja processado antes): $f"
        return 0
    fi

    # Ja esta no remoto com checksum batendo? (retomabilidade -- nao
    # resube a toa se uma rodada anterior ja deixou isso pronto)
    if rclone check "$f" "$REMOTE_DIR" --one-way --include "$(basename "$f")" >/dev/null 2>&1; then
        log "JA ESTAVA NO REMOTO (checksum confere) -- removendo do cluster sem re-upload: $f"
        rm -f "$f" && log "REMOVIDO do cluster: $f" || log "ERRO ao remover (mesmo ja verificado no remoto): $f"
        return 0
    fi

    log "Iniciando upload: $f ($(du -h "$f" | cut -f1))"

    if rclone copy "$f" "$REMOTE_DIR"; then
        log "Upload concluido, verificando integridade: $f"
    else
        log "FALHOU o upload (rclone copy retornou erro) -- MANTENDO arquivo local: $f"
        return 1
    fi

    # Verificacao real: compara checksum local x remoto, nao so "terminou
    # sem erro". So apaga se isso bater.
    if rclone check "$f" "$REMOTE_DIR" --one-way --include "$(basename "$f")" >/dev/null 2>&1; then
        rm -f "$f"
        log "VERIFICADO e REMOVIDO do cluster: $f"
    else
        log "AVISO: upload terminou mas checksum NAO confere -- MANTENDO arquivo local: $f"
        return 1
    fi
}

log "===== Iniciando lote de ${#FILES[@]} arquivo(s), $PARALLEL em paralelo ====="

# Processa em lotes de $PARALLEL, esperando cada lote terminar antes do
# proximo (evita abrir mais conexoes simultaneas do que o testado).
total=${#FILES[@]}
i=0
while [ $i -lt $total ]; do
    batch=("${FILES[@]:$i:$PARALLEL}")
    log "--- Lote: ${batch[*]}"
    for f in "${batch[@]}"; do
        upload_verify_delete "$f" &
    done
    wait
    i=$((i + PARALLEL))
done

log "===== Lote completo. Resumo final: ====="
for f in "${FILES[@]}"; do
    if [ -f "$f" ]; then
        log "  AINDA NO CLUSTER (verifique o motivo acima): $f"
    else
        log "  OK, removido: $f"
    fi
done
log "===== Fim ====="
