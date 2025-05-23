# Define constants for steps and L
$steps = 1024 * 32
$L = 64

# Define the specific N_species values to run
$N_species_values = @(5.01, 5.05, 5.06, 5.08, 5.1, 5.2, 5.3, 5.5, 5.9)

# Define the number of simulations per N value
$N_sims = 1

# Define the directory containing the executable
$scriptDir = "C:\GitHub\Sequential-Spatial-Killing\src\probabilityPuyoPuyo"

# Get the total number of logical processors and define the maximum number of concurrent jobs
$totalProcessors = (Get-WmiObject -Class Win32_ComputerSystem).NumberOfLogicalProcessors
$maxJobs = [math]::Max(1, $totalProcessors - 4)  # Leave 4 processors free

foreach ($N_species in $N_species_values) {
    for ($sim = 1; $sim -le $N_sims; $sim++) {
        while ((Get-Job | Where-Object { $_.State -eq 'Running' }).Count -ge $maxJobs) {
            Start-Sleep -Seconds 1
        }

        Start-Job -ScriptBlock {
            param($scriptDir, $L, $N_species, $steps)
            & "$scriptDir\onlyAvalanche2D.exe" $L $N_species $steps
            if ($LASTEXITCODE -ne 0) {
                Write-Error "Execution failed for N_species=$N_species"
            }
        } -ArgumentList $scriptDir, $L, $N_species, $steps
    }
}

Get-Job | Wait-Job
Get-Job | Remove-Job

Write-Host "All simulations completed."